# Overwrite policy version inherited vom ParaView CMake.
cmake_policy(PUSH)
cmake_policy(VERSION 3.12...3.31)

if (POLICY CMP0169)
  cmake_policy(SET CMP0169 OLD)
endif ()

find_package(CGAL CONFIG REQUIRED)

vtk_module_add_module(VofFlow::VofDataUtils FORCE_STATIC HEADER_DIRECTORIES
  NOWRAP_HEADERS
    Grid/DataInterpolation.h
    Grid/DomainInfo.h
    Grid/Gradient.h
    Grid/GridData.h
    Grid/GridDataSet.h
    Grid/GridIterator.h
    Grid/GridTypes.h
    Math/Vector.h
    Misc/CgalUtil.h
    Misc/CgalVariant.h
    Misc/ListSearch.h
    Misc/Profiling.h
    Misc/VofData.h
    Plic/CachedPlic.h
    Plic/Plic.h
    Plic/Plic3.h
    Plic/PlicPolyData.h
    Plic/PlicUtil.h
    Plic/PolyCell.h
    Plic/Polygon3D.h
    VtkUtil/VtkUtilArray.h
    VtkUtil/VtkUtilMPI.h
    VtkUtil/VtkUtilWriter.h
  SOURCES
    Grid/DataInterpolation.cpp
    Grid/DomainInfo.cpp
    Grid/Gradient.cpp
    Plic/Plic.cpp
    Plic/Plic3.cpp
    Plic/PlicPolyData.cpp
    Plic/PlicUtil.cpp
    Plic/PolyCell.cpp
    Plic/Polygon3D.cpp
    VtkUtil/VtkUtilArray.cpp
    VtkUtil/VtkUtilWriter.cpp)
target_compile_features(VofDataUtils PUBLIC cxx_std_17)
set_target_properties(VofDataUtils PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(VofDataUtils PUBLIC CGAL::CGAL)

option(VOFFLOW_USE_TRACY "" OFF)
if (VOFFLOW_USE_TRACY)
  # We require TRACY_DELAYED_INIT. This would require an vcpkg overlay port, use fetch content instead.
  FetchContent_Declare(tracy
    URL "https://github.com/wolfpld/tracy/archive/v0.11.1.tar.gz"
    URL_HASH SHA256=2c11ca816f2b756be2730f86b0092920419f3dabc7a7173829ffd897d91888a1)
  FetchContent_GetProperties(tracy)
  if (NOT tracy_POPULATED)
    message(STATUS "Fetch tracy ...")
    FetchContent_Populate(tracy)
    option(TRACY_ENABLE "" ON)
    option(TRACY_DELAYED_INIT "" OFF)
    option(TRACY_STATIC "" OFF)
    add_subdirectory(${tracy_SOURCE_DIR} ${tracy_BINARY_DIR}) # TODO EXCLUDE_FROM_ALL but need to install .so
    mark_as_advanced(FORCE
      FETCHCONTENT_SOURCE_DIR_TRACY
      FETCHCONTENT_UPDATES_DISCONNECTED_TRACY)
  endif ()

  target_link_libraries(VofDataUtils PUBLIC TracyClient)

  target_compile_definitions(VofDataUtils PUBLIC VOFFLOW_USE_TRACY)
endif ()

cmake_policy(POP)
