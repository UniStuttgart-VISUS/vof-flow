cmake_minimum_required(VERSION 3.12...3.29 FATAL_ERROR)

include(CMakeDependentOption)

option(VOFFLOW_USE_VCPKG "Use vcpkg to load dependencies." ON)
cmake_dependent_option(VOFFLOW_DISABLE_VCPKG_BOOST "Disable vcpkg boost for ParaView superbuild." OFF "VOFFLOW_USE_VCPKG" OFF)

# vcpkg
if (VOFFLOW_USE_VCPKG)
  include(FetchContent)
  mark_as_advanced(FORCE
    FETCHCONTENT_BASE_DIR
    FETCHCONTENT_FULLY_DISCONNECTED
    FETCHCONTENT_QUIET
    FETCHCONTENT_UPDATES_DISCONNECTED)

  # Require git for download
  find_package(Git REQUIRED)

  FetchContent_Declare(vcpkg-download
    GIT_REPOSITORY https://github.com/microsoft/vcpkg.git
    GIT_TAG 2024.04.26
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/vcpkg)
  FetchContent_GetProperties(vcpkg-download)
  if (NOT vcpkg-download_POPULATED)
    message(STATUS "Fetch vcpkg ...")
    FetchContent_Populate(vcpkg-download)
    mark_as_advanced(FORCE
      FETCHCONTENT_SOURCE_DIR_VCPKG-DOWNLOAD
      FETCHCONTENT_UPDATES_DISCONNECTED_VCPKG-DOWNLOAD)
  endif ()

  set(VCPKG_BOOTSTRAP_OPTIONS "-disableMetrics")
  set(VCPKG_INSTALL_OPTIONS "--clean-after-build" "--no-print-usage")
  set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_BINARY_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
  set(ENV{VCPKG_FORCE_DOWNLOADED_BINARIES} ON)

  if (WIN32)
    set(VCPKG_TARGET_TRIPLET "x64-windows-static")
  endif ()

  if (VOFFLOW_DISABLE_VCPKG_BOOST)
    include(cmake/disable_vcpkg_boost.cmake)
  endif ()
endif ()

project(VofFlow
  LANGUAGES C CXX)

# Set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif ()

find_package(ParaView REQUIRED)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

paraview_plugin_scan(
  PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/ParaViewPlugin/paraview.plugin"
  PROVIDES_PLUGINS plugins
  ENABLE_BY_DEFAULT ON)

paraview_plugin_build(
  RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY "${PARAVIEW_PLUGIN_SUBDIR}"
  INSTALL_HEADERS OFF
  PLUGINS ${plugins})

option(VOFFLOW_TOOLS "Build tools." ON)
if (VOFFLOW_TOOLS)
  add_subdirectory(src/Tools)
endif ()
