# vtkIsotropicRemeshingFilter

An isotropic surface remeshing filter for VTK.

It is an independent implementation of the incremental remeshing algorithm
published by Botsch & Kobbelt, *"A Remeshing Approach to Multiresolution
Modeling"* (Symposium on Geometry Processing 2004), built only on VTK data
structures and based solely on the published paper (no external geometry
processing libraries are used or referenced; BSD-license compatible).
Each iteration performs:

1. **Split** — split every edge longer than 4/3 of the target edge length at
   its midpoint
2. **Collapse** — collapse every edge shorter than 4/5 of the target edge
   length into its midpoint (suppressed when it would create an edge longer
   than 4/3 of the target, or break the mesh topology)
3. **Flip** — flip edges to drive vertex valences towards 6 (4 on the
   boundary)
4. **Tangential relaxation** — smooth by projecting the uniform Laplacian
   barycenter onto the vertex tangent plane
5. **Projection** — project vertices back onto the original surface with
   `vtkStaticCellLocator`

Internally a half-edge structure is built to carry out the local operations
(split / collapse / flip).

## Features and limitations

- **Boundary preservation**: boundary edges are always constrained. Boundary
  vertices only move along the boundary polyline, and geometric corners are
  kept fixed.
- **Feature protection** (optional): with `ProtectFeaturesOn()`, edges whose
  dihedral angle exceeds `FeatureAngle` are protected the same way as
  boundaries.
- The input must be a **manifold surface mesh with consistent orientation**.
  Non-triangular cells are triangulated and duplicate points are merged
  internally. For inconsistently oriented inputs, apply `vtkPolyDataNormals`
  (`ConsistencyOn`) beforehand.
- **Attribute interpolation** (`InterpolateAttributes`, on by default):
  point data is carried over by barycentric interpolation at the closest
  point on the original surface, and cell data is copied from the original
  cell closest to each output triangle centroid. Since the output vertices
  end up on the original surface after the final projection, the point data
  interpolation is effectively an exact linear interpolation.
- The relaxation and projection steps are parallelized with `vtkSMPTools`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `TargetEdgeLength` | `0` (auto) | Target edge length. If `<= 0`, the mean edge length of the input is used |
| `NumberOfIterations` | `3` | Number of remeshing iterations |
| `NumberOfRelaxationSteps` | `1` | Tangential relaxation steps per iteration |
| `DoProject` | `true` | Whether to project back onto the original surface |
| `InterpolateAttributes` | `true` | Whether to interpolate point/cell data onto the output |
| `ProtectFeatures` | `false` | Whether to protect sharp feature edges |
| `FeatureAngle` | `60` (degrees) | Dihedral angle threshold for feature/corner detection |

## Usage

```cpp
#include "vtkIsotropicRemeshingFilter.h"

vtkNew<vtkIsotropicRemeshingFilter> remesh;
remesh->SetInputData(surface);          // vtkPolyData
remesh->SetTargetEdgeLength(0.5);       // choose to match the model scale
remesh->SetNumberOfIterations(5);
remesh->ProtectFeaturesOn();            // to keep sharp edges of CAD-like models
remesh->SetFeatureAngle(60.0);
remesh->Update();
vtkPolyData* result = remesh->GetOutput();
```

## Building

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DVTK_DIR=<VTK build/install directory>
cmake --build build --config Release
ctest --test-dir build -C Release
```

## Command-line tool

`RemeshSurface` (built with `BUILD_EXAMPLES=ON`) remeshes a file directly.

```
RemeshSurface input.{vtp|stl|ply|obj} output.{vtp|stl|ply} [targetEdgeLength] [iterations] [featureAngle]
```

- When `targetEdgeLength` is omitted (or 0), the mean edge length of the
  input is used
- Passing `featureAngle` enables feature protection

## Directory layout

```
src/        vtkIsotropicRemeshingFilter (the library)
tests/      regression tests (sphere, bounded plane, cube feature protection),
            integrity-checking debug harness, per-phase benchmark
examples/   RemeshSurface command-line tool
```

## References

- M. Botsch, L. Kobbelt, "A Remeshing Approach to Multiresolution Modeling",
  Symposium on Geometry Processing, 2004.
