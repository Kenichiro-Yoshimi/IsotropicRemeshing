// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
/**
 * @class   vtkIsotropicRemeshingFilter
 * @brief   remesh a triangulated surface into a mesh with uniform edge lengths
 *
 * vtkIsotropicRemeshingFilter implements the incremental isotropic remeshing
 * algorithm published by Botsch & Kobbelt, "A Remeshing Approach to
 * Multiresolution Modeling" (Symposium on Geometry Processing 2004).
 * Each iteration performs:
 *
 * 1. Split all edges longer than 4/3 of the target edge length at their
 *    midpoint.
 * 2. Collapse all edges shorter than 4/5 of the target edge length into
 *    their midpoint (unless this would create an edge longer than 4/3 of
 *    the target length, or break the mesh topology).
 * 3. Flip edges to drive vertex valences towards 6 (4 on the boundary).
 * 4. Relax vertices tangentially (uniform Laplacian projected onto the
 *    vertex tangent plane).
 * 5. Project vertices back onto the original input surface.
 *
 * The input must be a manifold, consistently oriented surface mesh.
 * Non-triangular cells are triangulated and duplicate points are merged
 * internally before remeshing. Boundary edges are always preserved:
 * boundary vertices are only moved along the boundary polyline and
 * geometric boundary corners are kept fixed. When ProtectFeatures is
 * enabled, edges whose dihedral angle exceeds FeatureAngle are treated
 * the same way as boundary edges.
 *
 * When InterpolateAttributes is enabled (the default), point data is
 * interpolated onto the remeshed vertices with the barycentric coordinates
 * of their closest point on the original surface, and cell data is copied
 * from the original cell closest to each output triangle centroid.
 *
 * @warning The input must not contain non-manifold edges. Inconsistently
 * oriented inputs can be fixed beforehand with vtkPolyDataNormals
 * (ConsistencyOn).
 */

#ifndef vtkIsotropicRemeshingFilter_h
#define vtkIsotropicRemeshingFilter_h

#include "vtkPolyDataAlgorithm.h"

VTK_ABI_NAMESPACE_BEGIN

class vtkGenericCell;
class vtkStaticCellLocator;

class vtkIsotropicRemeshingFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkIsotropicRemeshingFilter* New();
  vtkTypeMacro(vtkIsotropicRemeshingFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  ///@{
  /**
   * Set/Get the target edge length of the remeshed output.
   * If the value is <= 0 (the default), the mean edge length of the
   * input mesh is used.
   */
  vtkSetMacro(TargetEdgeLength, double);
  vtkGetMacro(TargetEdgeLength, double);
  ///@}

  ///@{
  /**
   * Set/Get the number of remeshing iterations (split / collapse / flip /
   * relax / project rounds). More iterations give a more uniform mesh.
   * Values are clamped to be at least 1. Default is 3.
   */
  vtkSetClampMacro(NumberOfIterations, int, 1, VTK_INT_MAX);
  vtkGetMacro(NumberOfIterations, int);
  ///@}

  ///@{
  /**
   * Set/Get the number of tangential relaxation steps performed inside
   * each remeshing iteration. Values are clamped to be at least 0.
   * Default is 1.
   */
  vtkSetClampMacro(NumberOfRelaxationSteps, int, 0, VTK_INT_MAX);
  vtkGetMacro(NumberOfRelaxationSteps, int);
  ///@}

  ///@{
  /**
   * Enable/Disable projection of the vertices back onto the original
   * input surface after each iteration. Default is on.
   */
  vtkSetMacro(DoProject, bool);
  vtkGetMacro(DoProject, bool);
  vtkBooleanMacro(DoProject, bool);
  ///@}

  ///@{
  /**
   * Enable/Disable interpolation of the input attributes onto the remeshed
   * output. Point data is interpolated at the closest point on the original
   * surface; cell data is copied from the closest original cell.
   * Default is on.
   */
  vtkSetMacro(InterpolateAttributes, bool);
  vtkGetMacro(InterpolateAttributes, bool);
  vtkBooleanMacro(InterpolateAttributes, bool);
  ///@}

  ///@{
  /**
   * Enable/Disable protection of sharp feature edges. When enabled, edges
   * whose dihedral angle exceeds FeatureAngle are constrained: they are
   * never flipped, vertices on them only move along the feature polyline,
   * and feature corners are kept fixed. Default is off.
   */
  vtkSetMacro(ProtectFeatures, bool);
  vtkGetMacro(ProtectFeatures, bool);
  vtkBooleanMacro(ProtectFeatures, bool);
  ///@}

  ///@{
  /**
   * Set/Get the dihedral angle (in degrees) above which an edge is
   * considered a sharp feature when ProtectFeatures is enabled. The same
   * angle is used to detect corners on boundary / feature polylines.
   * Values are clamped to the range [0, 180]. Default is 60.
   */
  vtkSetClampMacro(FeatureAngle, double, 0.0, 180.0);
  vtkGetMacro(FeatureAngle, double);
  ///@}

protected:
  vtkIsotropicRemeshingFilter();
  ~vtkIsotropicRemeshingFilter() override = default;

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  /**
   * Interpolate the attributes of the original surface onto the remeshed
   * output (closest-point barycentric interpolation for point data, closest
   * cell for cell data).
   */
  void CopyAttributes(
    vtkPolyData* surface, vtkPolyData* output, vtkStaticCellLocator* locator, vtkGenericCell* cell);

  double TargetEdgeLength = 0.0;
  int NumberOfIterations = 3;
  int NumberOfRelaxationSteps = 1;
  bool DoProject = true;
  bool InterpolateAttributes = true;
  bool ProtectFeatures = false;
  double FeatureAngle = 60.0;

private:
  vtkIsotropicRemeshingFilter(const vtkIsotropicRemeshingFilter&) = delete;
  void operator=(const vtkIsotropicRemeshingFilter&) = delete;
};

VTK_ABI_NAMESPACE_END

#endif
