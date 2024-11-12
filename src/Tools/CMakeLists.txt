# Tools

find_package(OpenMP REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

add_executable(MaxVelocity MaxVelocity.cpp)
target_compile_features(MaxVelocity PUBLIC cxx_std_17)
set_target_properties(MaxVelocity PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(MaxVelocity
  PUBLIC
    ParaView::VTKExtensionsIOCore
    VofDataUtils
    OpenMP::OpenMP_CXX)
if (UNIX)
  set_target_properties(MaxVelocity PROPERTIES
    BUILD_RPATH_USE_ORIGIN TRUE
    INSTALL_RPATH "\$ORIGIN/../lib")
endif ()

add_executable(SeedGrid SeedGrid.cpp)
target_compile_features(SeedGrid PUBLIC cxx_std_17)
set_target_properties(SeedGrid PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(SeedGrid
  PUBLIC
    ParaView::VTKExtensionsIOCore
    nlohmann_json::nlohmann_json
    VofDataUtils
    VofTracking)
if (UNIX)
  set_target_properties(SeedGrid PROPERTIES
    BUILD_RPATH_USE_ORIGIN TRUE
    INSTALL_RPATH "\$ORIGIN/../lib")
endif ()

add_executable(PolyGap PolyGap.cpp)
target_compile_features(PolyGap PUBLIC cxx_std_17)
set_target_properties(PolyGap PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(PolyGap
  PUBLIC
    ParaView::VTKExtensionsIOCore
    VofDataUtils)
if (UNIX)
  set_target_properties(PolyGap PROPERTIES
    BUILD_RPATH_USE_ORIGIN TRUE
    INSTALL_RPATH "\$ORIGIN/../lib")
endif ()

add_executable(SmoothNormals SmoothNormals.cpp)
target_compile_features(SmoothNormals PUBLIC cxx_std_17)
set_target_properties(SmoothNormals PROPERTIES CXX_EXTENSIONS OFF)
target_link_libraries(SmoothNormals
  PUBLIC
    ParaView::VTKExtensionsIOCore
    VTK::FiltersParallel
    VofDataUtils)
if (UNIX)
  set_target_properties(SmoothNormals PROPERTIES
    BUILD_RPATH_USE_ORIGIN TRUE
    INSTALL_RPATH "\$ORIGIN/../lib")
endif ()

# Install
include(GNUInstallDirs)

install(TARGETS MaxVelocity SeedGrid PolyGap SmoothNormals
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
