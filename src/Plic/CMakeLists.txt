# Overwrite policy version inherited vom ParaView CMake.
cmake_policy(PUSH)
cmake_policy(VERSION 3.12...3.31)

vtk_module_add_module(VofFlow::Plic FORCE_STATIC
  CLASSES
    vtkAbstractPlic
    vtkPlic
    vtkPlic3)
target_compile_features(Plic PUBLIC cxx_std_17)
set_target_properties(Plic PROPERTIES CXX_EXTENSIONS OFF)

cmake_policy(POP)
