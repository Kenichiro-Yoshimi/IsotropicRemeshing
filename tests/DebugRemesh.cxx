// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
// Debug harness: runs the remeshing phases step by step with integrity
// checks after each phase. Includes the implementation file directly to
// reach the internal half-edge structures.

#include "../src/vtkIsotropicRemeshingFilter.cxx" // NOLINT

#include <vtkCleanPolyData.h>
#include <vtkNew.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>

#include <cstdlib>
#include <iostream>

namespace
{

bool Check(const HalfedgeMesh& mesh, const char* phase, int iter)
{
  std::string err = mesh.CheckIntegrity();
  if (!err.empty())
  {
    std::cerr << "INTEGRITY FAILURE after " << phase << " (iteration " << iter << "): " << err
              << std::endl;
    return false;
  }
  return true;
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

} // namespace

int main(int, char*[])
{
  vtkNew<vtkSphereSource> sphere;
  sphere->SetRadius(1.0);
  sphere->SetThetaResolution(60);
  sphere->SetPhiResolution(10);

  vtkNew<vtkCleanPolyData> clean;
  clean->SetInputConnection(sphere->GetOutputPort());
  vtkNew<vtkTriangleFilter> tri;
  tri->SetInputConnection(clean->GetOutputPort());
  tri->Update();

  HalfedgeMesh mesh;
  if (!mesh.Build(tri->GetOutput()))
  {
    std::cerr << "Build failed: " << mesh.ErrorMessage << std::endl;
    return EXIT_FAILURE;
  }
  mesh.DetectCorners(60.0);
  if (!Check(mesh, "Build", 0))
  {
    return EXIT_FAILURE;
  }
  std::cout << "built: " << LiveCount(mesh.VRemoved) << " vertices, " << LiveCount(mesh.FRemoved)
            << " faces" << std::endl;

  Remesher remesher(mesh, 0.1);
  for (int iter = 0; iter < 6; ++iter)
  {
    remesher.SplitLongEdges();
    std::cout << "iter " << iter << " after split:    V=" << LiveCount(mesh.VRemoved)
              << " F=" << LiveCount(mesh.FRemoved) << std::endl;
    if (!Check(mesh, "SplitLongEdges", iter))
    {
      return EXIT_FAILURE;
    }
    int collapsed = remesher.CollapseShortEdges();
    std::cout << "iter " << iter << " after collapse: V=" << LiveCount(mesh.VRemoved)
              << " F=" << LiveCount(mesh.FRemoved) << " (collapsed " << collapsed << ")"
              << std::endl;
    if (!Check(mesh, "CollapseShortEdges", iter))
    {
      return EXIT_FAILURE;
    }
    remesher.EqualizeValences();
    if (!Check(mesh, "EqualizeValences", iter))
    {
      return EXIT_FAILURE;
    }
    remesher.TangentialRelaxation();
    if (!Check(mesh, "TangentialRelaxation", iter))
    {
      return EXIT_FAILURE;
    }
  }
  std::cout << "all phases consistent" << std::endl;
  return EXIT_SUCCESS;
}
