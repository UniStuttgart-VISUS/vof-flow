# Disable all vcpkg boost* ports, to use the version provided by the ParaView superbuild.
# This overwrites all boost ports with empty ports. The empty ports use the original
# vcpkg.json to keep all features and dependencies next to an empty portfile.

set(empty_ports_dir "${CMAKE_CURRENT_BINARY_DIR}/vcpkg_empty_ports")

# Cleanup port dir
file(REMOVE_RECURSE "${empty_ports_dir}")
file(MAKE_DIRECTORY "${empty_ports_dir}")

# Searching boost ports
set(port_dirs "${CMAKE_CURRENT_BINARY_DIR}/vcpkg/ports")
file(GLOB boost_ports RELATIVE ${port_dirs} ${port_dirs}/boost*)

# Create empty port
foreach (portname ${boost_ports})
  file(COPY "${port_dirs}/${portname}/vcpkg.json" DESTINATION "${empty_ports_dir}/${portname}/")
  file(WRITE "${empty_ports_dir}/${portname}/portfile.cmake" "set(VCPKG_POLICY_EMPTY_PACKAGE enabled)\n")
endforeach ()

set(VCPKG_OVERLAY_PORTS "${empty_ports_dir};${VCPKG_OVERLAY_PORTS}")
