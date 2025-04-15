# Changelog

## Version 1.2.0

- Add support for non-uniform grids.

## Version 1.1.0

- Refactor component extraction to generate label names ordered by global grid position to improve reproducibility.
  Previously, the labels were ordered by MPI process ID and local grid position within a process, making the label
  ordering dependent on the number of MPI processes.
- Replace the SeedGrid, PolyGap, and SmoothNormals post-processing binary tools with a single VTK Filter to avoid
  running extra binaries. The new filter runs on the seed point output from the separation boundary filter.
- Replace the MaxVelocity binary tool with a VTK Filter to avoid running extra binaries. In addition, the new filter
  now also generates a velocity histogram and supports running with MPI.
- Add support for vtkImageData as input dataset.
- Add support for datasets with an extent not starting at (0, 0, 0).
- Fix a bug with gradient calculation for datasets with only a single cell in one dimension.
- Fix image data output for correct rendering within ParaView when using MPI.
- Require at least ParaView 5.12.
- Update docs and reproducibility scripts to ParaView 5.13.2.
- Add script for the post-processing filter to reproducibility docs.
- Fix plugin path in reproducibility scripts for Windows.
- Use VTK internal nlohmann json to avoid conflicts with ParaView 5.13.
- Add support for CGAL 6.x. Add a wrapper for the breaking changes to keep support for CGAL 5.x.
- Update vcpkg, tracy, and CMake version range.
- Cleanup and other minor fixes: Replace direct console output with VTK logging macros, refactor utility functions, fix
  CI script, fix clang-tidy suggestions, code formatting, ...

## Version 1.0.0

- Initial release.
