// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
// Regression tests for vtkIsotropicRemeshingFilter.
//
// Test 1: remesh an anisotropically tessellated sphere and verify that the
//         result is a closed manifold triangle mesh with near-uniform edge
//         lengths and preserved area.
// Test 2: remesh a planar rectangle and verify that the boundary polyline
//         (including its four corners) is preserved.

#include "vtkIsotropicRemeshingFilter.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFeatureEdges.h>
#include <vtkIntArray.h>
#include <vtkMassProperties.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkPlaneSource.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPolyDataWriter.h>

#include <string>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

namespace
{

int NumFailures = 0;

// When non-empty, input/output meshes are written there as .vtp files for
// visual inspection (e.g. in ParaView).
std::string OutputDir;

void WriteMesh(vtkPolyData* mesh, const std::string& name)
{
  if (OutputDir.empty())
  {
    return;
  }
  vtkNew<vtkXMLPolyDataWriter> writer;
  std::string path = OutputDir + "/" + name + ".vtp";
  writer->SetFileName(path.c_str());
  writer->SetInputData(mesh);
  writer->Write();
  std::cout << "  wrote " << path << std::endl;
}

#define TEST_ASSERT(cond, msg)                                                                     \
  do                                                                                               \
  {                                                                                                \
    if (!(cond))                                                                                   \
    {                                                                                              \
      std::cerr << "FAILED: " << msg << " (line " << __LINE__ << ")" << std::endl;                 \
      ++NumFailures;                                                                               \
    }                                                                                              \
  } while (false)

struct EdgeStats
{
  double Min = 1e300;
  double Max = 0;
  double Mean = 0;
  vtkIdType Count = 0;
};

EdgeStats ComputeEdgeStats(vtkPolyData* mesh)
{
  EdgeStats stats;
  std::map<std::pair<vtkIdType, vtkIdType>, bool> seen;
  vtkCellArray* polys = mesh->GetPolys();
  vtkIdType npts;
  const vtkIdType* pts;
  double total = 0;
  for (polys->InitTraversal(); polys->GetNextCell(npts, pts);)
  {
    for (vtkIdType i = 0; i < npts; ++i)
    {
      vtkIdType a = pts[i];
      vtkIdType b = pts[(i + 1) % npts];
      auto key = std::make_pair(std::min(a, b), std::max(a, b));
      if (!seen.emplace(key, true).second)
      {
        continue;
      }
      double pa[3], pb[3];
      mesh->GetPoint(a, pa);
      mesh->GetPoint(b, pb);
      double len = std::sqrt(vtkMath::Distance2BetweenPoints(pa, pb));
      stats.Min = std::min(stats.Min, len);
      stats.Max = std::max(stats.Max, len);
      total += len;
      ++stats.Count;
    }
  }
  stats.Mean = stats.Count ? total / stats.Count : 0;
  return stats;
}

vtkIdType CountEdges(vtkPolyData* mesh, bool boundary, bool nonManifold)
{
  vtkNew<vtkFeatureEdges> edges;
  edges->SetInputData(mesh);
  edges->BoundaryEdgesOff();
  edges->FeatureEdgesOff();
  edges->ManifoldEdgesOff();
  edges->NonManifoldEdgesOff();
  if (boundary)
  {
    edges->BoundaryEdgesOn();
  }
  if (nonManifold)
  {
    edges->NonManifoldEdgesOn();
  }
  edges->Update();
  return edges->GetOutput()->GetNumberOfCells();
}

//------------------------------------------------------------------------------
void TestSphere()
{
  std::cout << "== Sphere test ==" << std::endl;
  const double targetLength = 0.1;

  vtkNew<vtkSphereSource> sphere;
  sphere->SetRadius(1.0);
  sphere->SetThetaResolution(60);
  sphere->SetPhiResolution(10); // strongly anisotropic triangles
  sphere->Update();

  vtkNew<vtkMassProperties> massIn;
  massIn->SetInputConnection(sphere->GetOutputPort());
  massIn->Update();
  double areaIn = massIn->GetSurfaceArea();

  vtkNew<vtkIsotropicRemeshingFilter> remesh;
  remesh->SetInputConnection(sphere->GetOutputPort());
  remesh->SetTargetEdgeLength(targetLength);
  remesh->SetNumberOfIterations(6);
  remesh->Update();
  vtkPolyData* out = remesh->GetOutput();

  std::cout << "  input:  " << sphere->GetOutput()->GetNumberOfPoints() << " points, "
            << sphere->GetOutput()->GetNumberOfPolys() << " triangles" << std::endl;
  std::cout << "  output: " << out->GetNumberOfPoints() << " points, " << out->GetNumberOfPolys()
            << " triangles" << std::endl;
  WriteMesh(sphere->GetOutput(), "sphere_input");
  WriteMesh(out, "sphere_remeshed");

  TEST_ASSERT(out->GetNumberOfPolys() > 100, "output has a reasonable number of triangles");

  // all cells are triangles
  vtkCellArray* polys = out->GetPolys();
  vtkIdType npts;
  const vtkIdType* pts;
  bool allTris = true;
  for (polys->InitTraversal(); polys->GetNextCell(npts, pts);)
  {
    allTris = allTris && (npts == 3);
  }
  TEST_ASSERT(allTris, "all output cells are triangles");

  // closed manifold surface
  TEST_ASSERT(CountEdges(out, true, false) == 0, "output sphere is closed (no boundary edges)");
  TEST_ASSERT(CountEdges(out, false, true) == 0, "output sphere has no non-manifold edges");

  EdgeStats stats = ComputeEdgeStats(out);
  std::cout << "  edge lengths: min=" << stats.Min << " mean=" << stats.Mean << " max=" << stats.Max
            << " (target " << targetLength << ")" << std::endl;
  TEST_ASSERT(stats.Mean > 0.75 * targetLength && stats.Mean < 1.25 * targetLength,
    "mean edge length is within 25% of the target");
  TEST_ASSERT(stats.Max < 2.0 * targetLength, "no edge is longer than twice the target");
  TEST_ASSERT(stats.Min > 0.25 * targetLength, "no edge is shorter than a quarter of the target");

  // area preservation (projection keeps the shape on the sphere)
  vtkNew<vtkMassProperties> massOut;
  massOut->SetInputData(out);
  massOut->Update();
  double areaOut = massOut->GetSurfaceArea();
  std::cout << "  area: in=" << areaIn << " out=" << areaOut << std::endl;
  TEST_ASSERT(std::abs(areaOut - areaIn) / areaIn < 0.05, "surface area is preserved within 5%");
}

//------------------------------------------------------------------------------
void TestPlaneBoundary()
{
  std::cout << "== Plane boundary test ==" << std::endl;
  const double targetLength = 0.07;

  vtkNew<vtkPlaneSource> plane;
  plane->SetOrigin(0.0, 0.0, 0.0);
  plane->SetPoint1(1.0, 0.0, 0.0);
  plane->SetPoint2(0.0, 1.0, 0.0);
  plane->SetXResolution(25);
  plane->SetYResolution(3); // strongly anisotropic quads
  vtkNew<vtkTriangleFilter> tri;
  tri->SetInputConnection(plane->GetOutputPort());

  vtkNew<vtkIsotropicRemeshingFilter> remesh;
  remesh->SetInputConnection(tri->GetOutputPort());
  remesh->SetTargetEdgeLength(targetLength);
  remesh->SetNumberOfIterations(6);
  remesh->Update();
  vtkPolyData* out = remesh->GetOutput();

  std::cout << "  output: " << out->GetNumberOfPoints() << " points, " << out->GetNumberOfPolys()
            << " triangles" << std::endl;
  WriteMesh(tri->GetOutput(), "plane_input");
  WriteMesh(out, "plane_remeshed");

  TEST_ASSERT(CountEdges(out, false, true) == 0, "output plane has no non-manifold edges");

  EdgeStats stats = ComputeEdgeStats(out);
  std::cout << "  edge lengths: min=" << stats.Min << " mean=" << stats.Mean << " max=" << stats.Max
            << " (target " << targetLength << ")" << std::endl;
  TEST_ASSERT(stats.Mean > 0.7 * targetLength && stats.Mean < 1.3 * targetLength,
    "mean edge length is within 30% of the target");

  // every point must stay in the plane
  bool inPlane = true;
  for (vtkIdType i = 0; i < out->GetNumberOfPoints(); ++i)
  {
    double p[3];
    out->GetPoint(i, p);
    inPlane = inPlane && std::abs(p[2]) < 1e-9;
  }
  TEST_ASSERT(inPlane, "all output points stay in the z=0 plane");

  // boundary points must lie on the unit-square perimeter
  vtkNew<vtkFeatureEdges> edges;
  edges->SetInputData(out);
  edges->BoundaryEdgesOn();
  edges->FeatureEdgesOff();
  edges->ManifoldEdgesOff();
  edges->NonManifoldEdgesOff();
  edges->Update();
  vtkPolyData* boundary = edges->GetOutput();
  TEST_ASSERT(boundary->GetNumberOfCells() > 0, "plane output has a boundary");

  auto onPerimeter = [](const double p[3])
  {
    auto on = [](double v, double t) { return std::abs(v - t) < 1e-7; };
    bool inX = p[0] > -1e-7 && p[0] < 1 + 1e-7;
    bool inY = p[1] > -1e-7 && p[1] < 1 + 1e-7;
    return ((on(p[0], 0) || on(p[0], 1)) && inY) || ((on(p[1], 0) || on(p[1], 1)) && inX);
  };
  bool allOnPerimeter = true;
  for (vtkIdType i = 0; i < boundary->GetNumberOfPoints(); ++i)
  {
    double p[3];
    boundary->GetPoint(i, p);
    if (!onPerimeter(p))
    {
      std::cerr << "  off-perimeter boundary point: (" << p[0] << ", " << p[1] << ", " << p[2]
                << ")" << std::endl;
      allOnPerimeter = false;
    }
  }
  TEST_ASSERT(allOnPerimeter, "all boundary points lie on the square perimeter");

  // the four corners must be preserved exactly
  const double corners[4][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { 1, 1 } };
  for (auto& corner : corners)
  {
    bool found = false;
    for (vtkIdType i = 0; i < out->GetNumberOfPoints() && !found; ++i)
    {
      double p[3];
      out->GetPoint(i, p);
      found = std::abs(p[0] - corner[0]) < 1e-9 && std::abs(p[1] - corner[1]) < 1e-9;
    }
    TEST_ASSERT(found, "corner (" << corner[0] << ", " << corner[1] << ") is preserved");
  }
}

//------------------------------------------------------------------------------
void TestFeatureProtection()
{
  std::cout << "== Feature protection test (unit cube) ==" << std::endl;
  // Build a unit cube out of triangles (two per face).
  vtkNew<vtkPoints> points;
  const double v[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
    { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
  for (auto& p : v)
  {
    points->InsertNextPoint(p);
  }
  // outward-oriented quads, split into triangles
  const vtkIdType quads[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 },
    { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
  vtkNew<vtkCellArray> tris;
  for (auto& q : quads)
  {
    vtkIdType t0[3] = { q[0], q[1], q[2] };
    vtkIdType t1[3] = { q[0], q[2], q[3] };
    tris->InsertNextCell(3, t0);
    tris->InsertNextCell(3, t1);
  }
  vtkNew<vtkPolyData> cube;
  cube->SetPoints(points);
  cube->SetPolys(tris);

  const double targetLength = 0.2;
  vtkNew<vtkIsotropicRemeshingFilter> remesh;
  remesh->SetInputData(cube);
  remesh->SetTargetEdgeLength(targetLength);
  remesh->SetNumberOfIterations(8);
  remesh->ProtectFeaturesOn();
  remesh->SetFeatureAngle(60.0);
  remesh->Update();
  vtkPolyData* out = remesh->GetOutput();

  std::cout << "  output: " << out->GetNumberOfPoints() << " points, " << out->GetNumberOfPolys()
            << " triangles" << std::endl;
  WriteMesh(cube, "cube_input");
  WriteMesh(out, "cube_remeshed");

  TEST_ASSERT(CountEdges(out, true, false) == 0, "output cube is closed");
  TEST_ASSERT(CountEdges(out, false, true) == 0, "output cube has no non-manifold edges");
  TEST_ASSERT(out->GetNumberOfPolys() > 12, "the cube was refined");

  // all points must stay on the cube surface, i.e. at least one coordinate
  // is 0 or 1
  bool onSurface = true;
  for (vtkIdType i = 0; i < out->GetNumberOfPoints(); ++i)
  {
    double p[3];
    out->GetPoint(i, p);
    bool ok = false;
    for (int c = 0; c < 3; ++c)
    {
      ok = ok || std::abs(p[c]) < 1e-7 || std::abs(p[c] - 1) < 1e-7;
    }
    bool inside = true;
    for (int c = 0; c < 3; ++c)
    {
      inside = inside && p[c] > -1e-7 && p[c] < 1 + 1e-7;
    }
    if (!(ok && inside))
    {
      std::cerr << "  off-surface point: (" << p[0] << ", " << p[1] << ", " << p[2] << ")"
                << std::endl;
      onSurface = false;
    }
  }
  TEST_ASSERT(onSurface, "all points stay on the cube surface");

  // the 8 cube corners must be preserved
  for (auto& corner : v)
  {
    bool found = false;
    for (vtkIdType i = 0; i < out->GetNumberOfPoints() && !found; ++i)
    {
      double p[3];
      out->GetPoint(i, p);
      found = std::abs(p[0] - corner[0]) < 1e-9 && std::abs(p[1] - corner[1]) < 1e-9 &&
        std::abs(p[2] - corner[2]) < 1e-9;
    }
    TEST_ASSERT(found,
      "cube corner (" << corner[0] << ", " << corner[1] << ", " << corner[2] << ") is preserved");
  }
}

//------------------------------------------------------------------------------
void TestAttributeInterpolation()
{
  std::cout << "== Attribute interpolation test ==" << std::endl;
  vtkNew<vtkSphereSource> sphere;
  sphere->SetRadius(1.0);
  sphere->SetThetaResolution(60);
  sphere->SetPhiResolution(10);
  sphere->Update();
  vtkPolyData* input = sphere->GetOutput();

  // point scalar: the z coordinate (linear over each triangle, so it must be
  // reproduced almost exactly by the barycentric interpolation)
  vtkNew<vtkDoubleArray> elevation;
  elevation->SetName("Elevation");
  elevation->SetNumberOfTuples(input->GetNumberOfPoints());
  for (vtkIdType i = 0; i < input->GetNumberOfPoints(); ++i)
  {
    double p[3];
    input->GetPoint(i, p);
    elevation->SetValue(i, p[2]);
  }
  input->GetPointData()->AddArray(elevation);

  // cell scalar: hemisphere label by cell centroid
  vtkNew<vtkIntArray> hemisphere;
  hemisphere->SetName("Hemisphere");
  hemisphere->SetNumberOfTuples(input->GetNumberOfCells());
  for (vtkIdType i = 0; i < input->GetNumberOfCells(); ++i)
  {
    double bounds[6];
    input->GetCellBounds(i, bounds);
    double zMid = 0.5 * (bounds[4] + bounds[5]);
    hemisphere->SetValue(i, zMid > 0 ? 1 : 0);
  }
  input->GetCellData()->AddArray(hemisphere);

  vtkNew<vtkIsotropicRemeshingFilter> remesh;
  remesh->SetInputData(input);
  remesh->SetTargetEdgeLength(0.1);
  remesh->SetNumberOfIterations(6);
  remesh->Update();
  vtkPolyData* out = remesh->GetOutput();
  WriteMesh(input, "attributes_input");
  WriteMesh(out, "attributes_remeshed");

  auto* outElevation = vtkDoubleArray::SafeDownCast(out->GetPointData()->GetArray("Elevation"));
  TEST_ASSERT(outElevation != nullptr, "Elevation point array is present on the output");
  TEST_ASSERT(outElevation && outElevation->GetNumberOfTuples() == out->GetNumberOfPoints(),
    "Elevation array has one tuple per output point");
  if (outElevation)
  {
    double maxErr = 0;
    for (vtkIdType i = 0; i < out->GetNumberOfPoints(); ++i)
    {
      double p[3];
      out->GetPoint(i, p);
      maxErr = std::max(maxErr, std::abs(outElevation->GetValue(i) - p[2]));
    }
    std::cout << "  max |Elevation - z| = " << maxErr << std::endl;
    TEST_ASSERT(maxErr < 1e-6, "point data is interpolated at the projected positions");
  }

  auto* outHemisphere = vtkIntArray::SafeDownCast(out->GetCellData()->GetArray("Hemisphere"));
  TEST_ASSERT(outHemisphere != nullptr, "Hemisphere cell array is present on the output");
  TEST_ASSERT(outHemisphere && outHemisphere->GetNumberOfTuples() == out->GetNumberOfCells(),
    "Hemisphere array has one tuple per output cell");
  if (outHemisphere)
  {
    vtkIdType wrong = 0;
    for (vtkIdType i = 0; i < out->GetNumberOfCells(); ++i)
    {
      double bounds[6];
      out->GetCellBounds(i, bounds);
      double zMid = 0.5 * (bounds[4] + bounds[5]);
      // The input band straddling the equator (phi in [80, 100] degrees)
      // reaches |z| = cos(80 deg) ~ 0.174 and carries the "0" label, so an
      // exact label can only be expected safely above that band.
      if (std::abs(zMid) < 0.25)
      {
        continue;
      }
      if (outHemisphere->GetValue(i) != (zMid > 0 ? 1 : 0))
      {
        if (wrong < 5)
        {
          double b[6];
          out->GetCellBounds(i, b);
          std::cerr << "  cell " << i << " zMid=" << zMid << " label=" << outHemisphere->GetValue(i)
                    << " bounds z=[" << b[4] << "," << b[5] << "]" << std::endl;
        }
        ++wrong;
      }
    }
    std::cout << "  mislabeled cells away from the equator: " << wrong << std::endl;
    TEST_ASSERT(wrong == 0, "cell data is copied from the matching original cells");
  }
}

} // namespace

//------------------------------------------------------------------------------
// VTK-style test entry point, suitable for the vtkTestDriver when the class
// is moved into a VTK module.
int TestIsotropicRemeshingFilter(int argc, char* argv[])
{
  if (argc > 1)
  {
    OutputDir = argv[1];
  }
  TestSphere();
  TestPlaneBoundary();
  TestFeatureProtection();
  TestAttributeInterpolation();

  if (NumFailures)
  {
    std::cerr << NumFailures << " check(s) failed." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "All checks passed." << std::endl;
  return EXIT_SUCCESS;
}

//------------------------------------------------------------------------------
// Standalone entry point; a VTK module build generates this via the test
// driver instead.
int main(int argc, char* argv[])
{
  return TestIsotropicRemeshingFilter(argc, argv);
}
