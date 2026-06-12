cmake_minimum_required(VERSION 3.22)
project(IsotropicRemeshing CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(VTK 9 REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  CommonExecutionModel
  FiltersCore
)

add_library(IsotropicRemeshing STATIC
  src/vtkIsotropicRemeshingFilter.cxx
)
target_include_directories(IsotropicRemeshing PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
target_link_libraries(IsotropicRemeshing PUBLIC
  VTK::CommonCore
  VTK::CommonDataModel
  VTK::CommonExecutionModel
  VTK::FiltersCore
)

# ---------------------------------------------------------------------------
option(BUILD_TESTING "Build the tests" ON)
if(BUILD_TESTING)
  enable_testing()
  find_package(VTK 9 REQUIRED COMPONENTS FiltersSources IOXML)
  add_executable(TestIsotropicRemeshingFilter
    tests/TestIsotropicRemeshingFilter.cxx
  )
  target_link_libraries(TestIsotropicRemeshingFilter PRIVATE
    IsotropicRemeshing
    VTK::FiltersSources
    VTK::IOXML
  )
  add_test(NAME IsotropicRemeshingFilter COMMAND TestIsotropicRemeshingFilter)

  # step-by-step debugging harness with integrity checks
  add_executable(DebugRemesh tests/DebugRemesh.cxx)
  target_link_libraries(DebugRemesh PRIVATE
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::CommonExecutionModel
    VTK::FiltersCore
    VTK::FiltersSources
  )

  # per-phase timing benchmark
  add_executable(BenchmarkRemesh tests/BenchmarkRemesh.cxx)
  target_link_libraries(BenchmarkRemesh PRIVATE
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::CommonExecutionModel
    VTK::FiltersCore
    VTK::FiltersSources
  )
  if(WIN32)
    # let the test find the VTK DLLs
    set_tests_properties(IsotropicRemeshingFilter PROPERTIES
      ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${VTK_PREFIX_PATH}/bin/$<CONFIG>"
    )
  endif()
endif()

# ---------------------------------------------------------------------------
option(BUILD_EXAMPLES "Build the example command-line tool" ON)
if(BUILD_EXAMPLES)
  find_package(VTK 9 REQUIRED COMPONENTS IOGeometry IOPLY IOXML)
  add_executable(RemeshSurface examples/RemeshSurface.cxx)
  target_link_libraries(RemeshSurface PRIVATE
    IsotropicRemeshing
    VTK::IOGeometry
    VTK::IOPLY
    VTK::IOXML
  )
  vtk_module_autoinit(TARGETS RemeshSurface MODULES VTK::IOGeometry VTK::IOPLY VTK::IOXML)
endif()
