# Overwrite policy version inherited vom ParaView CMake.
cmake_policy(PUSH)
cmake_policy(VERSION 3.12...3.31)

find_package(OpenMP REQUIRED)
find_package(CGAL CONFIG REQUIRED)
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE CACHE INTERNAL "" FORCE)

vtk_module_add_module(VofFlow::VofTracking FORCE_STATIC
  CLASSES
    vtkVofTracking
    vtkVofBoundary
    vtkVofMaxVelocity
  SOURCES
    Advection.cpp
    Boundary.cpp
    ComponentExtraction.cpp
    Particles.cpp
    Seed.cpp
    Uncertainty.cpp)
target_compile_features(VofTracking PUBLIC cxx_std_17)
set_target_properties(VofTracking PROPERTIES CXX_EXTENSIONS OFF)

vtk_module_link(VofFlow::VofTracking
  PRIVATE
    CGAL::CGAL
    OpenMP::OpenMP_CXX)

cmake_policy(POP)
