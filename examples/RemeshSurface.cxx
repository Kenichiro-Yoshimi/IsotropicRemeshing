// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
// Command-line isotropic remeshing of a surface mesh file.
//
// Usage:
//   RemeshSurface input.{vtp|stl|ply|obj} output.{vtp|stl|ply}
//                 [targetEdgeLength] [iterations] [featureAngle]
//
// If targetEdgeLength is omitted or 0, the mean edge length of the input is
// used. If featureAngle is given, feature protection is enabled.

#include "vtkIsotropicRemeshingFilter.h"

#include <vtkNew.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

namespace
{

std::string Extension(const std::string& path)
{
  size_t dot = path.find_last_of('.');
  std::string ext = dot == std::string::npos ? "" : path.substr(dot + 1);
  std::transform(
    ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
  return ext;
}

vtkSmartPointer<vtkPolyData> ReadSurface(const std::string& path)
{
  std::string ext = Extension(path);
  vtkSmartPointer<vtkPolyData> surface;
  auto read = [&path, &surface](auto reader)
  {
    reader->SetFileName(path.c_str());
    reader->Update();
    surface = reader->GetOutput();
  };
  if (ext == "vtp")
  {
    read(vtkSmartPointer<vtkXMLPolyDataReader>::New());
  }
  else if (ext == "stl")
  {
    read(vtkSmartPointer<vtkSTLReader>::New());
  }
  else if (ext == "ply")
  {
    read(vtkSmartPointer<vtkPLYReader>::New());
  }
  else if (ext == "obj")
  {
    read(vtkSmartPointer<vtkOBJReader>::New());
  }
  return surface;
}

bool WriteSurface(vtkPolyData* surface, const std::string& path)
{
  std::string ext = Extension(path);
  auto write = [&path, surface](auto writer)
  {
    writer->SetFileName(path.c_str());
    writer->SetInputData(surface);
    return writer->Write() == 1;
  };
  if (ext == "vtp")
  {
    return write(vtkSmartPointer<vtkXMLPolyDataWriter>::New());
  }
  if (ext == "stl")
  {
    return write(vtkSmartPointer<vtkSTLWriter>::New());
  }
  if (ext == "ply")
  {
    return write(vtkSmartPointer<vtkPLYWriter>::New());
  }
  return false;
}

} // namespace

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0]
              << " input.{vtp|stl|ply|obj} output.{vtp|stl|ply}"
                 " [targetEdgeLength] [iterations] [featureAngle]"
              << std::endl;
    return EXIT_FAILURE;
  }
  double targetLength = argc > 3 ? std::atof(argv[3]) : 0.0;
  int iterations = argc > 4 ? std::atoi(argv[4]) : 5;
  double featureAngle = argc > 5 ? std::atof(argv[5]) : -1.0;

  vtkSmartPointer<vtkPolyData> surface = ReadSurface(argv[1]);
  if (!surface || surface->GetNumberOfPoints() == 0)
  {
    std::cerr << "Could not read a surface from " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Input:  " << surface->GetNumberOfPoints() << " points, "
            << surface->GetNumberOfCells() << " cells" << std::endl;

  vtkNew<vtkIsotropicRemeshingFilter> remesh;
  remesh->SetInputData(surface);
  remesh->SetTargetEdgeLength(targetLength);
  remesh->SetNumberOfIterations(iterations);
  if (featureAngle >= 0)
  {
    remesh->ProtectFeaturesOn();
    remesh->SetFeatureAngle(featureAngle);
  }
  remesh->Update();
  vtkPolyData* output = remesh->GetOutput();
  if (output->GetNumberOfPolys() == 0)
  {
    std::cerr << "Remeshing failed." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Output: " << output->GetNumberOfPoints() << " points, "
            << output->GetNumberOfPolys() << " triangles" << std::endl;

  if (!WriteSurface(output, argv[2]))
  {
    std::cerr << "Could not write " << argv[2] << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Wrote " << argv[2] << std::endl;
  return EXIT_SUCCESS;
}
