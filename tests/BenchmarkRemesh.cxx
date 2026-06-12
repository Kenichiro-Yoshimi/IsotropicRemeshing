// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
// Benchmark harness: times each remeshing phase on inputs of increasing size.
// Includes the implementation file directly to reach the internal phases.

#include "../src/vtkIsotropicRemeshingFilter.cxx" // NOLINT

#include <vtkCleanPolyData.h>
#include <vtkNew.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>

#include <chrono>
#include <cstdlib>
#include <iostream>

namespace
{

using Clock = std::chrono::steady_clock;

double Ms(Clock::time_point t0, Clock::time_point t1)
{
  return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

size_t LiveCount(const std::vector<bool>& removed)
{
  size_t n = 0;
  for (bool r : removed)
  {
    n += r ? 0 : 1;
  }
  return n;
}

void RunCase(int thetaRes, int phiRes, double targetLength, int iterations)
{
  vtkNew<vtkSphereSource> sphere;
  sphere->SetRadius(1.0);
  sphere->SetThetaResolution(thetaRes);
  sphere->SetPhiResolution(phiRes);
  vtkNew<vtkCleanPolyData> clean;
  clean->SetInputConnection(sphere->GetOutputPort());
  vtkNew<vtkTriangleFilter> tri;
  tri->SetInputConnection(clean->GetOutputPort());
  tri->Update();
  vtkPolyData* surface = tri->GetOutput();

  std::cout << "=== input " << surface->GetNumberOfPoints() << " pts / "
            << surface->GetNumberOfPolys() << " tris, target " << targetLength << ", " << iterations
            << " iterations ===" << std::endl;

  auto tStart = Clock::now();
  HalfedgeMesh mesh;
  if (!mesh.Build(surface))
  {
    std::cerr << "Build failed" << std::endl;
    return;
  }
  mesh.DetectCorners(60.0);
  auto tBuild = Clock::now();
  std::cout << "  build halfedge:  " << Ms(tStart, tBuild) << " ms" << std::endl;

  vtkNew<vtkStaticCellLocator> locator;
  locator->SetDataSet(surface);
  locator->BuildLocator();
  auto tLoc = Clock::now();
  std::cout << "  build locator:   " << Ms(tBuild, tLoc) << " ms" << std::endl;

  double split = 0, collapse = 0, flip = 0, relax = 0, project = 0;
  Remesher remesher(mesh, targetLength);
  for (int iter = 0; iter < iterations; ++iter)
  {
    auto t0 = Clock::now();
    remesher.SplitLongEdges();
    auto t1 = Clock::now();
    remesher.CollapseShortEdges();
    auto t2 = Clock::now();
    remesher.EqualizeValences();
    auto t3 = Clock::now();
    remesher.TangentialRelaxation();
    auto t4 = Clock::now();
    remesher.ProjectToSurface(locator);
    auto t5 = Clock::now();
    split += Ms(t0, t1);
    collapse += Ms(t1, t2);
    flip += Ms(t2, t3);
    relax += Ms(t3, t4);
    project += Ms(t4, t5);
  }
  auto tEnd = Clock::now();
  std::cout << "  split:           " << split << " ms" << std::endl;
  std::cout << "  collapse:        " << collapse << " ms" << std::endl;
  std::cout << "  flip:            " << flip << " ms" << std::endl;
  std::cout << "  relax:           " << relax << " ms" << std::endl;
  std::cout << "  project:         " << project << " ms" << std::endl;
  std::cout << "  TOTAL:           " << Ms(tStart, tEnd) << " ms, output "
            << LiveCount(mesh.VRemoved) << " pts / " << LiveCount(mesh.FRemoved) << " tris"
            << std::endl;
}

} // namespace

int main(int, char*[])
{
  RunCase(60, 30, 0.05, 5);    //  ~3.5k tris in,  ~12k tris out
  RunCase(200, 100, 0.02, 5);  //  ~40k tris in,   ~73k tris out
  RunCase(400, 200, 0.01, 5);  // ~160k tris in,  ~290k tris out
  RunCase(800, 400, 0.005, 5); // ~640k tris in, ~1.16M tris out
  return EXIT_SUCCESS;
}
