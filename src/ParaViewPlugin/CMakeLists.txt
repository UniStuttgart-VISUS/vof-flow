# Overwrite policy version inherited vom ParaView CMake.
cmake_policy(PUSH)
cmake_policy(VERSION 3.12...3.29)

paraview_add_plugin(VofFlow
  VERSION "1.0"
  SERVER_MANAGER_XML VofFlow.xml
  MODULES
    VofFlow::VofDataUtils
    VofFlow::VofTracking
    VofFlow::Plic
  MODULE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/../VofDataUtils/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/../VofTracking/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/../Plic/vtk.module"
  UI_RESOURCES VofFlow.qrc)

target_link_libraries(VofFlow
  PRIVATE
    VofFlow::VofDataUtils
    VofFlow::VofTracking
    VofFlow::Plic)

cmake_policy(POP)
