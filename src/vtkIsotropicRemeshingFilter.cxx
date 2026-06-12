// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkIsotropicRemeshingFilter.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSMPThreadLocalObject.h"
#include "vtkSMPTools.h"
#include "vtkStaticCellLocator.h"
#include "vtkTriangleFilter.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{

using Vec3 = std::array<double, 3>;

inline Vec3 operator-(const Vec3& a, const Vec3& b)
{
  return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}
inline Vec3 operator+(const Vec3& a, const Vec3& b)
{
  return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}
inline Vec3 operator*(double s, const Vec3& a)
{
  return { s * a[0], s * a[1], s * a[2] };
}
inline double Dot(const Vec3& a, const Vec3& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vec3 Cross(const Vec3& a, const Vec3& b)
{
  return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}
inline double Norm(const Vec3& a)
{
  return std::sqrt(Dot(a, a));
}

//------------------------------------------------------------------------------
// A compact half-edge mesh supporting the local operations required by
// incremental isotropic remeshing: edge split, edge collapse and edge flip.
// Boundaries are represented by explicit half-edges with Face == -1, linked
// into boundary loops, so Next/Prev are defined for every half-edge.
struct HalfedgeMesh
{
  // Vertices
  std::vector<Vec3> Points;
  std::vector<int> VHalfedge; // an outgoing half-edge (boundary one if any)
  std::vector<bool> VRemoved;
  std::vector<bool> VCorner; // corner vertices are never moved or removed

  // Half-edges
  std::vector<int> HVertex; // vertex the half-edge points to
  std::vector<int> HNext;
  std::vector<int> HPrev;
  std::vector<int> HOpp;
  std::vector<int> HFace; // -1 for boundary half-edges
  std::vector<bool> HRemoved;
  std::vector<bool> HConstrained; // boundary or protected feature edge

  // Faces
  std::vector<int> FHalfedge;
  std::vector<bool> FRemoved;

  std::string ErrorMessage;

  // Set when a circulation hits its iteration cap, which indicates a
  // corrupted connectivity; checked by CheckIntegrity(). Atomic because
  // circulations may run concurrently during relaxation and projection.
  mutable std::atomic<bool> CirculationOverflow{ false };

  int Source(int h) const { return this->HVertex[this->HOpp[h]]; }
  int Target(int h) const { return this->HVertex[h]; }
  bool IsBoundaryEdge(int h) const { return this->HFace[h] < 0 || this->HFace[this->HOpp[h]] < 0; }

  double EdgeLength(int h) const
  {
    return Norm(this->Points[this->Target(h)] - this->Points[this->Source(h)]);
  }

  Vec3 FaceNormal(int f) const
  {
    int h = this->FHalfedge[f];
    const Vec3& a = this->Points[this->Source(h)];
    const Vec3& b = this->Points[this->Target(h)];
    const Vec3& c = this->Points[this->Target(this->HNext[h])];
    return Cross(b - a, c - a);
  }

  // Visit the outgoing half-edges of vertex v. The iteration count is capped
  // so that a corrupted mesh cannot hang the circulation; the overflow flag
  // is checked by CheckIntegrity().
  template <typename Functor>
  void CirculateOutgoing(int v, Functor&& f) const
  {
    int start = this->VHalfedge[v];
    if (start < 0)
    {
      return;
    }
    int h = start;
    size_t guard = 0;
    do
    {
      if (!f(h))
      {
        return;
      }
      h = this->HNext[this->HOpp[h]];
      if (++guard > this->HVertex.size())
      {
        this->CirculationOverflow = true;
        return;
      }
    } while (h != start);
  }

  bool IsBoundaryVertex(int v) const
  {
    bool boundary = false;
    this->CirculateOutgoing(v,
      [this, &boundary](int h)
      {
        if (this->HFace[h] < 0 || this->HFace[this->HOpp[h]] < 0)
        {
          boundary = true;
          return false;
        }
        return true;
      });
    return boundary;
  }

  int Valence(int v) const
  {
    int n = 0;
    this->CirculateOutgoing(v,
      [&n](int)
      {
        ++n;
        return true;
      });
    return n;
  }

  // Number of constrained edges incident to v.
  int ConstrainedDegree(int v) const
  {
    int n = 0;
    this->CirculateOutgoing(v,
      [this, &n](int h)
      {
        if (this->HConstrained[h])
        {
          ++n;
        }
        return true;
      });
    return n;
  }

  // A vertex lying on a constrained polyline (boundary or feature).
  bool IsConstrainedVertex(int v) const { return this->ConstrainedDegree(v) > 0; }

  bool AreNeighbors(int v, int w) const
  {
    bool found = false;
    this->CirculateOutgoing(v,
      [this, &found, w](int h)
      {
        if (this->Target(h) == w)
        {
          found = true;
          return false;
        }
        return true;
      });
    return found;
  }

  // Make VHalfedge[v] point to a boundary half-edge if the vertex is on the
  // boundary, so that circulation visits the full one-ring.
  void AdjustOutgoing(int v)
  {
    int start = this->VHalfedge[v];
    if (start < 0)
    {
      return;
    }
    int h = start;
    do
    {
      if (this->HFace[h] < 0)
      {
        this->VHalfedge[v] = h;
        return;
      }
      h = this->HNext[this->HOpp[h]];
    } while (h != start);
  }

  //----------------------------------------------------------------------------
  // Validate the half-edge connectivity. Returns an empty string when the
  // structure is consistent. Intended for tests and debugging.
  std::string CheckIntegrity() const
  {
    this->CirculationOverflow = false;
    auto fail = [](const std::string& msg) { return msg; };
    int numHe = static_cast<int>(this->HVertex.size());
    for (int h = 0; h < numHe; ++h)
    {
      if (this->HRemoved[h])
      {
        continue;
      }
      int o = this->HOpp[h];
      if (o < 0 || o >= numHe || o == h || this->HRemoved[o] || this->HOpp[o] != h)
      {
        return fail("bad opposite at halfedge " + std::to_string(h));
      }
      if (this->HNext[h] < 0 || this->HRemoved[this->HNext[h]] ||
        this->HPrev[this->HNext[h]] != h || this->HRemoved[this->HPrev[h]] ||
        this->HNext[this->HPrev[h]] != h)
      {
        return fail("bad next/prev at halfedge " + std::to_string(h));
      }
      int v = this->HVertex[h];
      if (v < 0 || this->VRemoved[v])
      {
        return fail("halfedge " + std::to_string(h) + " points to a removed vertex");
      }
      if (this->HVertex[o] == v)
      {
        return fail("self-loop edge at halfedge " + std::to_string(h));
      }
      if (this->HConstrained[h] != this->HConstrained[o])
      {
        return fail("constraint mismatch at halfedge " + std::to_string(h));
      }
      int f = this->HFace[h];
      if (f >= 0 && this->FRemoved[f])
      {
        return fail("halfedge " + std::to_string(h) + " belongs to a removed face");
      }
      if (this->HNext[this->HNext[this->HNext[h]]] != h && f >= 0)
      {
        return fail("face loop of halfedge " + std::to_string(h) + " is not a triangle");
      }
      if (f >= 0 && (this->HFace[this->HNext[h]] != f || this->HFace[this->HPrev[h]] != f))
      {
        return fail("inconsistent face at halfedge " + std::to_string(h));
      }
    }
    for (size_t f = 0; f < this->FHalfedge.size(); ++f)
    {
      if (!this->FRemoved[f] &&
        (this->FHalfedge[f] < 0 || this->HRemoved[this->FHalfedge[f]] ||
          this->HFace[this->FHalfedge[f]] != static_cast<int>(f)))
      {
        return fail("bad face -> halfedge link at face " + std::to_string(f));
      }
    }
    for (size_t v = 0; v < this->Points.size(); ++v)
    {
      if (this->VRemoved[v] || this->VHalfedge[v] < 0)
      {
        continue;
      }
      int h = this->VHalfedge[v];
      if (this->HRemoved[h])
      {
        return fail("vertex " + std::to_string(v) + " points to a removed halfedge");
      }
      if (this->Source(h) != static_cast<int>(v))
      {
        return fail("outgoing halfedge of vertex " + std::to_string(v) + " has wrong source");
      }
      int valence = this->Valence(static_cast<int>(v));
      if (this->CirculationOverflow)
      {
        return fail("circulation around vertex " + std::to_string(v) + " does not terminate");
      }
      if (valence < 2)
      {
        return fail("vertex " + std::to_string(v) + " has valence " + std::to_string(valence));
      }
    }
    return {};
  }

  //----------------------------------------------------------------------------
  bool Build(vtkPolyData* mesh)
  {
    vtkIdType numPts = mesh->GetNumberOfPoints();
    vtkIdType numTris = mesh->GetNumberOfPolys();
    this->Points.resize(numPts);
    for (vtkIdType i = 0; i < numPts; ++i)
    {
      double p[3];
      mesh->GetPoint(i, p);
      this->Points[i] = { p[0], p[1], p[2] };
    }
    this->VHalfedge.assign(numPts, -1);
    this->VRemoved.assign(numPts, false);
    this->VCorner.assign(numPts, false);

    this->HVertex.reserve(6 * numTris);
    this->FHalfedge.reserve(numTris);

    // Map a directed edge (a -> b) to its half-edge, to detect opposites and
    // non-manifold configurations.
    std::unordered_map<long long, int> directed;
    directed.reserve(3 * numTris);
    auto key = [numPts](int a, int b) { return static_cast<long long>(a) * numPts + b; };

    vtkCellArray* polys = mesh->GetPolys();
    vtkIdType npts;
    const vtkIdType* pts;
    for (polys->InitTraversal(); polys->GetNextCell(npts, pts);)
    {
      if (npts != 3)
      {
        this->ErrorMessage = "Internal triangulation failed: non-triangular cell encountered.";
        return false;
      }
      int a = static_cast<int>(pts[0]);
      int b = static_cast<int>(pts[1]);
      int c = static_cast<int>(pts[2]);
      if (a == b || b == c || c == a)
      {
        continue; // skip degenerate triangles
      }
      int f = static_cast<int>(this->FHalfedge.size());
      int h0 = static_cast<int>(this->HVertex.size());
      const int verts[3] = { a, b, c };
      for (int i = 0; i < 3; ++i)
      {
        int src = verts[i];
        int dst = verts[(i + 1) % 3];
        if (!directed.emplace(key(src, dst), h0 + i).second)
        {
          this->ErrorMessage = "The input mesh is non-manifold or inconsistently oriented "
                               "(a directed edge is used by more than one triangle). "
                               "vtkPolyDataNormals with ConsistencyOn may fix the orientation.";
          return false;
        }
        this->HVertex.push_back(dst);
        this->HNext.push_back(h0 + (i + 1) % 3);
        this->HPrev.push_back(h0 + (i + 2) % 3);
        this->HOpp.push_back(-1);
        this->HFace.push_back(f);
        this->VHalfedge[src] = h0 + i;
      }
      this->FHalfedge.push_back(h0);
    }
    this->HRemoved.assign(this->HVertex.size(), false);
    this->HConstrained.assign(this->HVertex.size(), false);
    this->FRemoved.assign(this->FHalfedge.size(), false);

    // Link opposite half-edges; create boundary half-edges for unmatched ones.
    int numInterior = static_cast<int>(this->HVertex.size());
    std::unordered_map<int, int> boundaryOut; // source vertex -> boundary half-edge
    for (int h = 0; h < numInterior; ++h)
    {
      if (this->HOpp[h] >= 0)
      {
        continue;
      }
      int a = this->HVertex[this->HPrev[h]]; // source of h (HOpp is not linked yet)
      int b = this->HVertex[h];
      auto it = directed.find(key(b, a));
      if (it != directed.end())
      {
        this->HOpp[h] = it->second;
        this->HOpp[it->second] = h;
      }
      else
      {
        // boundary half-edge b -> a
        int bh = static_cast<int>(this->HVertex.size());
        this->HVertex.push_back(a);
        this->HNext.push_back(-1);
        this->HPrev.push_back(-1);
        this->HOpp.push_back(h);
        this->HFace.push_back(-1);
        this->HRemoved.push_back(false);
        this->HConstrained.push_back(false);
        this->HOpp[h] = bh;
        if (!boundaryOut.emplace(b, bh).second)
        {
          this->ErrorMessage =
            "The input mesh is non-manifold (a vertex has more than one boundary fan).";
          return false;
        }
      }
    }
    // Link boundary loops: Next of a boundary half-edge is the boundary
    // half-edge starting at its target vertex.
    for (auto& kv : boundaryOut)
    {
      int bh = kv.second;
      auto it = boundaryOut.find(this->HVertex[bh]);
      if (it == boundaryOut.end())
      {
        this->ErrorMessage = "The input mesh has an inconsistent boundary (non-manifold).";
        return false;
      }
      this->HNext[bh] = it->second;
      this->HPrev[it->second] = bh;
      // prefer boundary half-edges as the stored outgoing half-edge
      this->VHalfedge[kv.first] = bh;
    }
    // Boundary edges are always constrained.
    for (int h = 0; h < static_cast<int>(this->HVertex.size()); ++h)
    {
      if (this->HFace[h] < 0)
      {
        this->HConstrained[h] = true;
        this->HConstrained[this->HOpp[h]] = true;
      }
    }
    return true;
  }

  //----------------------------------------------------------------------------
  // Mark edges whose dihedral angle exceeds featureAngle (degrees) as
  // constrained.
  void DetectFeatureEdges(double featureAngle)
  {
    double cosAngle = std::cos(vtkMath::RadiansFromDegrees(featureAngle));
    int numHe = static_cast<int>(this->HVertex.size());
    for (int h = 0; h < numHe; ++h)
    {
      int o = this->HOpp[h];
      if (h > o || this->HRemoved[h] || this->HConstrained[h] || this->HFace[h] < 0 ||
        this->HFace[o] < 0)
      {
        continue;
      }
      Vec3 n0 = this->FaceNormal(this->HFace[h]);
      Vec3 n1 = this->FaceNormal(this->HFace[o]);
      double l0 = Norm(n0);
      double l1 = Norm(n1);
      if (l0 < 1e-300 || l1 < 1e-300)
      {
        continue;
      }
      if (Dot(n0, n1) / (l0 * l1) < cosAngle)
      {
        this->HConstrained[h] = true;
        this->HConstrained[o] = true;
      }
    }
  }

  //----------------------------------------------------------------------------
  // Mark vertices that must never move: vertices with a constrained degree
  // other than 2, and vertices where the two incident constrained edges meet
  // at a sharp turn.
  void DetectCorners(double featureAngle)
  {
    double cosTurn = std::cos(vtkMath::RadiansFromDegrees(featureAngle));
    int numV = static_cast<int>(this->Points.size());
    for (int v = 0; v < numV; ++v)
    {
      if (this->VRemoved[v] || this->VHalfedge[v] < 0)
      {
        continue;
      }
      std::vector<int> constrainedNbrs;
      this->CirculateOutgoing(v,
        [this, &constrainedNbrs](int h)
        {
          if (this->HConstrained[h])
          {
            constrainedNbrs.push_back(this->Target(h));
          }
          return true;
        });
      if (constrainedNbrs.empty())
      {
        continue;
      }
      if (constrainedNbrs.size() != 2)
      {
        this->VCorner[v] = true;
        continue;
      }
      Vec3 d0 = this->Points[constrainedNbrs[0]] - this->Points[v];
      Vec3 d1 = this->Points[constrainedNbrs[1]] - this->Points[v];
      double l0 = Norm(d0);
      double l1 = Norm(d1);
      if (l0 < 1e-300 || l1 < 1e-300)
      {
        continue;
      }
      // straight polyline: d0 and d1 point in opposite directions
      if (Dot(d0, d1) / (l0 * l1) > -cosTurn)
      {
        this->VCorner[v] = true;
      }
    }
  }

  //----------------------------------------------------------------------------
  // Split the edge of h at position p. Returns the index of the new vertex.
  int SplitEdge(int h, const Vec3& p)
  {
    int o = this->HOpp[h];
    int a = this->Source(h);
    bool constrained = this->HConstrained[h];

    int m = static_cast<int>(this->Points.size());
    this->Points.push_back(p);
    this->VHalfedge.push_back(-1);
    this->VRemoved.push_back(false);
    this->VCorner.push_back(false);

    // Insert the new edge a-m so that h becomes m->b and o becomes b->m.
    int e1 = this->NewHalfedge();  // a -> m
    int e1o = this->NewHalfedge(); // m -> a
    this->HVertex[e1] = m;
    this->HVertex[e1o] = a;
    this->HOpp[e1] = e1o;
    this->HOpp[e1o] = e1;
    this->HConstrained[e1] = constrained;
    this->HConstrained[e1o] = constrained;

    int hp = this->HPrev[h];
    int onx = this->HNext[o];
    this->Link(hp, e1);
    this->Link(e1, h);
    this->HFace[e1] = this->HFace[h];
    this->Link(o, e1o);
    this->Link(e1o, onx);
    this->HFace[e1o] = this->HFace[o];
    this->HVertex[o] = m;

    this->VHalfedge[m] = h; // m -> b
    if (this->VHalfedge[a] == h)
    {
      this->VHalfedge[a] = e1;
    }

    // Triangulate the two incident faces (when present) by connecting m to
    // the opposite apexes.
    if (this->HFace[h] >= 0)
    {
      // quad loop: e1 (a->m), h (m->b), hn (b->c), hp (c->a)
      int hn = this->HNext[h];
      int f0 = this->HFace[h];
      int d1 = this->NewHalfedge();  // m -> c
      int d1o = this->NewHalfedge(); // c -> m
      this->HVertex[d1] = this->Target(hn);
      this->HVertex[d1o] = m;
      this->HOpp[d1] = d1o;
      this->HOpp[d1o] = d1;
      int f2 = this->NewFace();
      // face f2: e1 (a->m), d1 (m->c), hp (c->a)
      this->Link(e1, d1);
      this->Link(d1, hp);
      this->Link(hp, e1);
      this->HFace[e1] = f2;
      this->HFace[d1] = f2;
      this->HFace[hp] = f2;
      this->FHalfedge[f2] = d1;
      // face f0: h (m->b), hn (b->c), d1o (c->m)
      this->Link(hn, d1o);
      this->Link(d1o, h);
      this->HFace[d1o] = f0;
      this->FHalfedge[f0] = h;
    }
    if (this->HFace[o] >= 0)
    {
      // quad loop: o (b->m), e1o (m->a), on (a->d), op (d->b)
      int on = this->HNext[e1o];
      int op = this->HPrev[o];
      int f1 = this->HFace[o];
      int d2 = this->NewHalfedge();  // m -> d
      int d2o = this->NewHalfedge(); // d -> m
      this->HVertex[d2] = this->Target(on);
      this->HVertex[d2o] = m;
      this->HOpp[d2] = d2o;
      this->HOpp[d2o] = d2;
      int f3 = this->NewFace();
      // face f1: o (b->m), d2 (m->d), op (d->b)
      this->Link(o, d2);
      this->Link(d2, op);
      this->Link(op, o);
      this->HFace[d2] = f1;
      this->FHalfedge[f1] = o;
      // face f3: e1o (m->a), on (a->d), d2o (d->m)
      this->Link(on, d2o);
      this->Link(d2o, e1o);
      this->HFace[e1o] = f3;
      this->HFace[on] = f3;
      this->HFace[d2o] = f3;
      this->FHalfedge[f3] = d2o;
    }
    this->AdjustOutgoing(m);
    this->AdjustOutgoing(a);
    return m;
  }

  //----------------------------------------------------------------------------
  // Check whether collapsing h (removing its source vertex into its target
  // vertex) keeps the mesh manifold.
  bool IsCollapseOk(int h) const
  {
    int o = this->HOpp[h];
    int v0 = this->Source(h);
    int v1 = this->Target(h);

    int c = this->HFace[h] >= 0 ? this->Target(this->HNext[h]) : -1;
    int d = this->HFace[o] >= 0 ? this->Target(this->HNext[o]) : -1;
    if (c >= 0 && c == d)
    {
      return false;
    }
    // Collapsing a triangle whose two other edges are both on the boundary
    // would create a dangling edge.
    if (this->HFace[h] >= 0 && this->IsBoundaryEdge(this->HNext[h]) &&
      this->IsBoundaryEdge(this->HPrev[h]))
    {
      return false;
    }
    if (this->HFace[o] >= 0 && this->IsBoundaryEdge(this->HNext[o]) &&
      this->IsBoundaryEdge(this->HPrev[o]))
    {
      return false;
    }
    // An edge between two boundary vertices must itself be a boundary edge.
    if (this->IsBoundaryVertex(v0) && this->IsBoundaryVertex(v1) && !this->IsBoundaryEdge(h))
    {
      return false;
    }
    // Link condition: the one-rings of v0 and v1 may only intersect in the
    // apex vertices c and d. The rings are small (~6 vertices), so a linear
    // scan over a flat array beats a hash set here.
    std::vector<int> ring1;
    ring1.reserve(16);
    this->CirculateOutgoing(v1,
      [this, &ring1](int hh)
      {
        ring1.push_back(this->Target(hh));
        return true;
      });
    bool ok = true;
    this->CirculateOutgoing(v0,
      [this, &ring1, &ok, v1, c, d](int hh)
      {
        int w = this->Target(hh);
        if (w != v1 && w != c && w != d && std::find(ring1.begin(), ring1.end(), w) != ring1.end())
        {
          ok = false;
          return false;
        }
        return true;
      });
    return ok;
  }

  //----------------------------------------------------------------------------
  // Collapse the edge of h: the source vertex of h is removed, all its edges
  // are reconnected to the target vertex. IsCollapseOk(h) must hold.
  void Collapse(int h)
  {
    int o = this->HOpp[h];
    int v0 = this->Source(h);
    int v1 = this->Target(h);
    int hn = this->HNext[h];
    int hp = this->HPrev[h];
    int on = this->HNext[o];
    int op = this->HPrev[o];
    int fh = this->HFace[h];
    int fo = this->HFace[o];

    // Reroute all half-edges pointing to v0 so that they point to v1.
    std::vector<int> incoming;
    this->CirculateOutgoing(v0,
      [this, &incoming](int hh)
      {
        incoming.push_back(this->HOpp[hh]);
        return true;
      });
    for (int hi : incoming)
    {
      this->HVertex[hi] = v1;
    }

    // Unlink h and o from their loops.
    this->Link(hp, hn);
    this->Link(op, on);
    if (fh >= 0)
    {
      this->FHalfedge[fh] = hn;
    }
    if (fo >= 0)
    {
      this->FHalfedge[fo] = on;
    }
    this->VHalfedge[v1] = hn; // hn starts at v1
    this->HRemoved[h] = true;
    this->HRemoved[o] = true;
    this->VRemoved[v0] = true;
    this->VHalfedge[v0] = -1;

    // The loops of the (former) incident triangles now have length two;
    // remove them by gluing their outer half-edges together.
    if (fh >= 0)
    {
      this->CollapseLoop(hn);
    }
    if (fo >= 0)
    {
      this->CollapseLoop(on);
    }
    this->AdjustOutgoing(v1);
  }

  //----------------------------------------------------------------------------
  // Flip the (interior) edge of h. Validity must be checked by the caller.
  void Flip(int h)
  {
    int o = this->HOpp[h];
    int hn = this->HNext[h];
    int hp = this->HPrev[h];
    int on = this->HNext[o];
    int op = this->HPrev[o];
    int f0 = this->HFace[h];
    int f1 = this->HFace[o];
    int a = this->Source(h);
    int b = this->Target(h);
    int c = this->Target(hn);
    int d = this->Target(on);

    this->HVertex[h] = c; // h becomes d -> c
    this->HVertex[o] = d; // o becomes c -> d

    // face f0: on (a->d), h (d->c), hp (c->a)
    this->Link(on, h);
    this->Link(h, hp);
    this->Link(hp, on);
    this->HFace[on] = f0;
    this->FHalfedge[f0] = h;
    // face f1: op (d->b), hn (b->c), o (c->d)
    this->Link(op, hn);
    this->Link(hn, o);
    this->Link(o, op);
    this->HFace[hn] = f1;
    this->FHalfedge[f1] = o;

    if (this->VHalfedge[a] == h)
    {
      this->VHalfedge[a] = on;
    }
    if (this->VHalfedge[b] == o)
    {
      this->VHalfedge[b] = hn;
    }
  }

  //----------------------------------------------------------------------------
  bool IsFlipOk(int h) const
  {
    int o = this->HOpp[h];
    if (this->HRemoved[h] || this->HFace[h] < 0 || this->HFace[o] < 0 || this->HConstrained[h])
    {
      return false;
    }
    int c = this->Target(this->HNext[h]);
    int d = this->Target(this->HNext[o]);
    if (c == d || this->AreNeighbors(c, d))
    {
      return false;
    }
    // Reject flips that would create degenerate or folded triangles.
    int a = this->Source(h);
    int b = this->Target(h);
    const Vec3& pa = this->Points[a];
    const Vec3& pb = this->Points[b];
    const Vec3& pc = this->Points[c];
    const Vec3& pd = this->Points[d];
    Vec3 nOld = Cross(pb - pa, pc - pa) + Cross(pa - pb, pd - pb);
    Vec3 n0 = Cross(pd - pa, pc - pa); // triangle (a, d, c)
    Vec3 n1 = Cross(pb - pd, pc - pd); // triangle (d, b, c)
    double scale = Norm(nOld);
    if (scale < 1e-300 || Dot(n0, nOld) <= 1e-12 * scale * Norm(n0) ||
      Dot(n1, nOld) <= 1e-12 * scale * Norm(n1))
    {
      return false;
    }
    return true;
  }

  //----------------------------------------------------------------------------
  void ToPolyData(vtkPolyData* output) const
  {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    std::vector<vtkIdType> remap(this->Points.size(), -1);
    for (size_t v = 0; v < this->Points.size(); ++v)
    {
      if (!this->VRemoved[v] && this->VHalfedge[v] >= 0)
      {
        remap[v] = points->InsertNextPoint(this->Points[v].data());
      }
    }
    vtkNew<vtkCellArray> tris;
    for (size_t f = 0; f < this->FHalfedge.size(); ++f)
    {
      if (this->FRemoved[f])
      {
        continue;
      }
      int h = this->FHalfedge[f];
      vtkIdType ids[3] = { remap[this->Source(h)], remap[this->Target(h)],
        remap[this->Target(this->HNext[h])] };
      tris->InsertNextCell(3, ids);
    }
    output->SetPoints(points);
    output->SetPolys(tris);
  }

private:
  int NewHalfedge()
  {
    int h = static_cast<int>(this->HVertex.size());
    this->HVertex.push_back(-1);
    this->HNext.push_back(-1);
    this->HPrev.push_back(-1);
    this->HOpp.push_back(-1);
    this->HFace.push_back(-1);
    this->HRemoved.push_back(false);
    this->HConstrained.push_back(false);
    return h;
  }

  int NewFace()
  {
    int f = static_cast<int>(this->FHalfedge.size());
    this->FHalfedge.push_back(-1);
    this->FRemoved.push_back(false);
    return f;
  }

  void Link(int h0, int h1)
  {
    this->HNext[h0] = h1;
    this->HPrev[h1] = h0;
  }

  // Remove a face whose loop has degenerated to two half-edges (h and
  // Next[h]) by gluing their outer opposites together.
  void CollapseLoop(int h)
  {
    int h1 = this->HNext[h];
    int o = this->HOpp[h];
    int o1 = this->HOpp[h1];
    int va = this->HVertex[h];
    int vb = this->HVertex[h1];
    int f = this->HFace[h];

    this->HOpp[o] = o1;
    this->HOpp[o1] = o;
    bool constrained = this->HConstrained[h] || this->HConstrained[h1];
    this->HConstrained[o] = this->HConstrained[o] || constrained;
    this->HConstrained[o1] = this->HConstrained[o1] || constrained;

    if (f >= 0)
    {
      this->FRemoved[f] = true;
    }
    this->HRemoved[h] = true;
    this->HRemoved[h1] = true;

    if (this->VHalfedge[va] == h1)
    {
      this->VHalfedge[va] = o;
    }
    if (this->VHalfedge[vb] == h)
    {
      this->VHalfedge[vb] = o1;
    }
    this->AdjustOutgoing(va);
    this->AdjustOutgoing(vb);
  }
};

//------------------------------------------------------------------------------
// The remeshing driver implementing Botsch & Kobbelt 2004.
struct Remesher
{
  HalfedgeMesh& Mesh;
  double High; // 4/3 * target length
  double Low;  // 4/5 * target length

  Remesher(HalfedgeMesh& mesh, double targetLength)
    : Mesh(mesh)
    , High(4.0 / 3.0 * targetLength)
    , Low(4.0 / 5.0 * targetLength)
  {
  }

  //----------------------------------------------------------------------------
  void SplitLongEdges()
  {
    HalfedgeMesh& m = this->Mesh;
    std::vector<int> stack;
    int numHe = static_cast<int>(m.HVertex.size());
    for (int h = 0; h < numHe; ++h)
    {
      if (!m.HRemoved[h] && h < m.HOpp[h] && m.EdgeLength(h) > this->High)
      {
        stack.push_back(h);
      }
    }
    while (!stack.empty())
    {
      int h = stack.back();
      stack.pop_back();
      if (m.HRemoved[h] || m.EdgeLength(h) <= this->High)
      {
        continue;
      }
      int a = m.Source(h);
      int b = m.Target(h);
      Vec3 mid = 0.5 * (m.Points[a] + m.Points[b]);
      int v = m.SplitEdge(h, mid);
      // Only the two halves of the split edge are re-enqueued: their length
      // strictly halves, which guarantees termination. The diagonals created
      // by the face splits may still be long (e.g. in fans of skinny
      // triangles, where they do not shrink under repeated splitting); they
      // are dealt with in the next remeshing iteration, after flips and
      // relaxation have improved the local triangle shapes.
      m.CirculateOutgoing(v,
        [this, &m, &stack, a, b](int hh)
        {
          int t = m.Target(hh);
          if ((t == a || t == b) && m.EdgeLength(hh) > this->High)
          {
            stack.push_back(hh);
          }
          return true;
        });
    }
  }

  //----------------------------------------------------------------------------
  // Returns the number of collapsed edges.
  int CollapseShortEdges()
  {
    HalfedgeMesh& m = this->Mesh;
    std::vector<int> candidates;
    int numHe = static_cast<int>(m.HVertex.size());
    for (int h = 0; h < numHe; ++h)
    {
      if (!m.HRemoved[h] && h < m.HOpp[h] && m.EdgeLength(h) < this->Low)
      {
        candidates.push_back(h);
      }
    }
    int collapsed = 0;
    for (int h : candidates)
    {
      if (m.HRemoved[h] || m.EdgeLength(h) >= this->Low)
      {
        continue;
      }
      if (this->TryCollapse(h) || this->TryCollapse(m.HOpp[h]))
      {
        ++collapsed;
      }
    }
    return collapsed;
  }

  //----------------------------------------------------------------------------
  // Try to collapse h (removing its source vertex). Returns true on success.
  bool TryCollapse(int h)
  {
    HalfedgeMesh& m = this->Mesh;
    int v0 = m.Source(h);
    int v1 = m.Target(h);
    if (m.VCorner[v0])
    {
      return false; // corners are never removed
    }
    bool c0 = m.IsConstrainedVertex(v0);
    bool c1 = m.IsConstrainedVertex(v1);
    bool edgeConstrained = m.HConstrained[h];

    Vec3 target;
    if (c0 && !edgeConstrained)
    {
      // moving a constrained vertex off its polyline is not allowed
      return false;
    }
    if (c0 && edgeConstrained)
    {
      // both endpoints on the same constrained polyline
      target = m.VCorner[v1] ? m.Points[v1] : 0.5 * (m.Points[v0] + m.Points[v1]);
    }
    else if (c1)
    {
      target = m.Points[v1];
    }
    else
    {
      target = 0.5 * (m.Points[v0] + m.Points[v1]);
    }

    // Do not create edges longer than the split threshold.
    if (this->CreatesLongEdge(v0, v0, v1, target) || this->CreatesLongEdge(v1, v0, v1, target) ||
      !m.IsCollapseOk(h))
    {
      return false;
    }
    // Reject collapses that fold over surrounding triangles.
    if (this->CollapseCreatesFoldover(h, target))
    {
      return false;
    }
    m.Points[v1] = target;
    m.Collapse(h);
    return true;
  }

  //----------------------------------------------------------------------------
  // Check whether moving the collapse survivor to `target` would create an
  // edge longer than the split threshold among the neighbors of v.
  bool CreatesLongEdge(int v, int v0, int v1, const Vec3& target) const
  {
    HalfedgeMesh& m = this->Mesh;
    bool tooLong = false;
    m.CirculateOutgoing(v,
      [this, &m, &tooLong, &target, v0, v1](int hh)
      {
        int w = m.Target(hh);
        if (w != v0 && w != v1 && Norm(m.Points[w] - target) > this->High)
        {
          tooLong = true;
          return false;
        }
        return true;
      });
    return tooLong;
  }

  //----------------------------------------------------------------------------
  bool CollapseCreatesFoldover(int h, const Vec3& target) const
  {
    HalfedgeMesh& m = this->Mesh;
    int v0 = m.Source(h);
    int v1 = m.Target(h);
    int f0 = m.HFace[h];
    int f1 = m.HFace[m.HOpp[h]];
    return this->FoldoverAround(v0, v0, v1, f0, f1, target) ||
      this->FoldoverAround(v1, v0, v1, f0, f1, target);
  }

  //----------------------------------------------------------------------------
  // Check whether moving vertex v to `target` inverts any of its incident
  // faces, ignoring the two faces (f0, f1) that the collapse removes.
  bool FoldoverAround(int v, int v0, int v1, int f0, int f1, const Vec3& target) const
  {
    HalfedgeMesh& m = this->Mesh;
    const Vec3& saved = m.Points[v];
    bool foldover = false;
    m.CirculateOutgoing(v,
      [&m, &foldover, &saved, &target, v0, v1, f0, f1](int hh)
      {
        int f = m.HFace[hh];
        if (f < 0 || f == f0 || f == f1)
        {
          return true;
        }
        int a = m.Target(hh);
        int b = m.Target(m.HNext[hh]);
        if (a == v0 || a == v1 || b == v0 || b == v1)
        {
          return true; // face vanishes or is checked via the other vertex
        }
        Vec3 before = Cross(m.Points[a] - saved, m.Points[b] - saved);
        Vec3 after = Cross(m.Points[a] - target, m.Points[b] - target);
        if (Dot(before, after) <= 0)
        {
          foldover = true;
          return false;
        }
        return true;
      });
    return foldover;
  }

  //----------------------------------------------------------------------------
  void EqualizeValences()
  {
    HalfedgeMesh& m = this->Mesh;
    int numV = static_cast<int>(m.Points.size());
    int numHe = static_cast<int>(m.HVertex.size());
    // Cache valences and boundary flags once per sweep instead of circulating
    // around four vertices for every edge; flips update the cache locally.
    std::vector<int> valence(numV, 0);
    std::vector<bool> boundary(numV, false);
    for (int h = 0; h < numHe; ++h)
    {
      if (m.HRemoved[h])
      {
        continue;
      }
      ++valence[m.HVertex[h]];
      if (m.HFace[h] < 0)
      {
        boundary[m.HVertex[h]] = true;
        boundary[m.HVertex[m.HOpp[h]]] = true;
      }
    }
    for (int h = 0; h < numHe; ++h)
    {
      if (m.HRemoved[h] || h > m.HOpp[h] || !m.IsFlipOk(h))
      {
        continue;
      }
      int o = m.HOpp[h];
      int a = m.Source(h);
      int b = m.Target(h);
      int c = m.Target(m.HNext[h]);
      int d = m.Target(m.HNext[o]);
      auto deviation = [&valence, &boundary](int v, int delta)
      {
        int target = boundary[v] ? 4 : 6;
        int dev = valence[v] + delta - target;
        return dev * dev;
      };
      long before = deviation(a, 0) + deviation(b, 0) + deviation(c, 0) + deviation(d, 0);
      long after = deviation(a, -1) + deviation(b, -1) + deviation(c, 1) + deviation(d, 1);
      if (after < before)
      {
        m.Flip(h);
        --valence[a];
        --valence[b];
        ++valence[c];
        ++valence[d];
      }
    }
  }

  //----------------------------------------------------------------------------
  void TangentialRelaxation()
  {
    HalfedgeMesh& m = this->Mesh;
    int numV = static_cast<int>(m.Points.size());
    std::vector<Vec3> newPos(m.Points);
    // RelaxVertex only reads the mesh and writes newPos[v], so the vertices
    // can be processed concurrently.
    vtkSMPTools::For(0, numV,
      [this, &newPos](vtkIdType begin, vtkIdType end)
      {
        for (vtkIdType v = begin; v < end; ++v)
        {
          this->RelaxVertex(static_cast<int>(v), newPos);
        }
      });
    m.Points.swap(newPos);
  }

  //----------------------------------------------------------------------------
  // Move vertex v to the projection of its uniform Laplacian barycenter onto
  // its tangent plane; the result is written to newPos[v].
  void RelaxVertex(int v, std::vector<Vec3>& newPos) const
  {
    HalfedgeMesh& m = this->Mesh;
    if (m.VRemoved[v] || m.VHalfedge[v] < 0 || m.VCorner[v] || m.IsConstrainedVertex(v))
    {
      return;
    }
    // uniform Laplacian barycenter and area-weighted vertex normal
    Vec3 q = { 0, 0, 0 };
    Vec3 n = { 0, 0, 0 };
    int count = 0;
    m.CirculateOutgoing(v,
      [&m, &q, &n, &count](int hh)
      {
        q = q + m.Points[m.Target(hh)];
        if (m.HFace[hh] >= 0)
        {
          n = n + m.FaceNormal(m.HFace[hh]);
        }
        ++count;
        return true;
      });
    if (count == 0)
    {
      return;
    }
    q = (1.0 / count) * q;
    double len = Norm(n);
    if (len < 1e-300)
    {
      return;
    }
    n = (1.0 / len) * n;
    newPos[v] = q + Dot(n, m.Points[v] - q) * n;
  }

  //----------------------------------------------------------------------------
  void ProjectToSurface(vtkStaticCellLocator* locator)
  {
    HalfedgeMesh& m = this->Mesh;
    int numV = static_cast<int>(m.Points.size());
    // vtkStaticCellLocator queries are thread-safe after BuildLocator(); each
    // thread needs its own scratch cell. ProjectVertex writes only its own
    // Points[v].
    vtkSMPThreadLocalObject<vtkGenericCell> threadCell;
    vtkSMPTools::For(0, numV,
      [this, locator, &threadCell](vtkIdType begin, vtkIdType end)
      {
        vtkGenericCell* cell = threadCell.Local();
        for (vtkIdType v = begin; v < end; ++v)
        {
          this->ProjectVertex(static_cast<int>(v), locator, cell);
        }
      });
  }

  //----------------------------------------------------------------------------
  // Move vertex v to its closest point on the original surface.
  void ProjectVertex(int v, vtkStaticCellLocator* locator, vtkGenericCell* cell) const
  {
    HalfedgeMesh& m = this->Mesh;
    if (m.VRemoved[v] || m.VHalfedge[v] < 0 || m.VCorner[v] || m.IsConstrainedVertex(v))
    {
      return;
    }
    double closest[3];
    vtkIdType cellId;
    int subId;
    double dist2;
    locator->FindClosestPoint(m.Points[v].data(), closest, cell, cellId, subId, dist2);
    if (cellId >= 0)
    {
      m.Points[v] = { closest[0], closest[1], closest[2] };
    }
  }
};

} // anonymous namespace

VTK_ABI_NAMESPACE_BEGIN

vtkStandardNewMacro(vtkIsotropicRemeshingFilter);

//------------------------------------------------------------------------------
vtkIsotropicRemeshingFilter::vtkIsotropicRemeshingFilter() = default;

//------------------------------------------------------------------------------
int vtkIsotropicRemeshingFilter::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkPolyData* input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);

  if (!input || input->GetNumberOfPoints() == 0)
  {
    vtkErrorMacro("Empty input.");
    return 0;
  }

  // Merge duplicate points and triangulate.
  vtkNew<vtkCleanPolyData> clean;
  clean->SetInputData(input);
  clean->PointMergingOn();
  vtkNew<vtkTriangleFilter> triangulate;
  triangulate->SetInputConnection(clean->GetOutputPort());
  triangulate->PassLinesOff();
  triangulate->PassVertsOff();
  triangulate->Update();
  vtkPolyData* surface = triangulate->GetOutput();

  if (surface->GetNumberOfPolys() == 0)
  {
    vtkErrorMacro("The input does not contain any polygons.");
    return 0;
  }

  HalfedgeMesh mesh;
  if (!mesh.Build(surface))
  {
    vtkErrorMacro(<< mesh.ErrorMessage);
    return 0;
  }
  if (this->ProtectFeatures)
  {
    mesh.DetectFeatureEdges(this->FeatureAngle);
  }
  mesh.DetectCorners(this->FeatureAngle);

  // Determine the target edge length.
  double targetLength = this->TargetEdgeLength;
  if (targetLength <= 0)
  {
    double total = 0;
    vtkIdType count = 0;
    for (int h = 0; h < static_cast<int>(mesh.HVertex.size()); ++h)
    {
      if (h < mesh.HOpp[h])
      {
        total += mesh.EdgeLength(h);
        ++count;
      }
    }
    if (count == 0)
    {
      vtkErrorMacro("The input does not contain any edges.");
      return 0;
    }
    targetLength = total / count;
    vtkDebugMacro("Using mean input edge length as target edge length: " << targetLength);
  }

  // Locator on the original surface, used for the projection step and for
  // attribute interpolation.
  vtkNew<vtkStaticCellLocator> locator;
  vtkNew<vtkGenericCell> cell;
  if (this->DoProject || this->InterpolateAttributes)
  {
    locator->SetDataSet(surface);
    locator->BuildLocator();
  }

  Remesher remesher(mesh, targetLength);
  bool aborted = false;
  for (int iter = 0; iter < this->NumberOfIterations; ++iter)
  {
    if (this->CheckAbort())
    {
      aborted = true;
      break;
    }
    remesher.SplitLongEdges();
    remesher.CollapseShortEdges();
    remesher.EqualizeValences();
    for (int step = 0; step < this->NumberOfRelaxationSteps; ++step)
    {
      remesher.TangentialRelaxation();
    }
    if (this->DoProject)
    {
      remesher.ProjectToSurface(locator);
    }
    this->UpdateProgress((iter + 1.0) / this->NumberOfIterations);
  }
  if (!aborted)
  {
    mesh.ToPolyData(output);
    if (this->InterpolateAttributes)
    {
      this->CopyAttributes(surface, output, locator, cell);
    }
  }
  return 1;
}

//------------------------------------------------------------------------------
void vtkIsotropicRemeshingFilter::CopyAttributes(
  vtkPolyData* surface, vtkPolyData* output, vtkStaticCellLocator* locator, vtkGenericCell* cell)
{
  vtkPointData* inPD = surface->GetPointData();
  vtkCellData* inCD = surface->GetCellData();
  vtkPointData* outPD = output->GetPointData();
  vtkCellData* outCD = output->GetCellData();

  // Point data: interpolate with the barycentric coordinates of the closest
  // point on the original surface. After the final projection step the
  // output vertices lie on that surface, so this amounts to an exact linear
  // interpolation at the vertex position.
  vtkIdType numPts = output->GetNumberOfPoints();
  outPD->InterpolateAllocate(inPD, numPts);
  double closest[3], pcoords[3], weights[3];
  vtkIdType cellId;
  int subId;
  double dist2;
  for (vtkIdType i = 0; i < numPts; ++i)
  {
    double x[3];
    output->GetPoint(i, x);
    locator->FindClosestPoint(x, closest, cell, cellId, subId, dist2);
    if (cellId < 0 || cell->GetNumberOfPoints() != 3)
    {
      continue;
    }
    cell->EvaluatePosition(closest, nullptr, subId, pcoords, dist2, weights);
    outPD->InterpolatePoint(inPD, i, cell->GetPointIds(), weights);
  }

  // Cell data: copy from the original cell closest to the triangle centroid.
  vtkCellArray* polys = output->GetPolys();
  outCD->CopyAllocate(inCD, polys->GetNumberOfCells());
  vtkIdType npts;
  const vtkIdType* pts;
  vtkIdType outCellId = 0;
  for (polys->InitTraversal(); polys->GetNextCell(npts, pts); ++outCellId)
  {
    double centroid[3] = { 0, 0, 0 };
    for (vtkIdType j = 0; j < npts; ++j)
    {
      double p[3];
      output->GetPoint(pts[j], p);
      centroid[0] += p[0] / npts;
      centroid[1] += p[1] / npts;
      centroid[2] += p[2] / npts;
    }
    locator->FindClosestPoint(centroid, closest, cell, cellId, subId, dist2);
    if (cellId >= 0)
    {
      outCD->CopyData(inCD, cellId, outCellId);
    }
  }
}

//------------------------------------------------------------------------------
void vtkIsotropicRemeshingFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "TargetEdgeLength: " << this->TargetEdgeLength << "\n";
  os << indent << "NumberOfIterations: " << this->NumberOfIterations << "\n";
  os << indent << "NumberOfRelaxationSteps: " << this->NumberOfRelaxationSteps << "\n";
  os << indent << "DoProject: " << (this->DoProject ? "On" : "Off") << "\n";
  os << indent << "InterpolateAttributes: " << (this->InterpolateAttributes ? "On" : "Off") << "\n";
  os << indent << "ProtectFeatures: " << (this->ProtectFeatures ? "On" : "Off") << "\n";
  os << indent << "FeatureAngle: " << this->FeatureAngle << "\n";
}
VTK_ABI_NAMESPACE_END
