#!/bin/bash

# This script is based on the ParaViewEasyPluginBuilder
# https://gitlab.kitware.com/paraview/paraview-easy-plugin-builder

# Create output dir if not exists
if [[ -e plugin_build ]]; then
  echo "Error: Directory 'plugin_build' already exists. Please remove!"
  exit 1
fi
mkdir plugin_build

# Create build script running in the docker container
cat > "plugin_build/run.sh" <<EOF
#!/bin/bash

# Create dirs
mkdir -p /vofflow/src
mkdir -p /vofflow/build

# Download repo
git clone https://github.com/UniStuttgart-VISUS/vof-flow.git /vofflow/src

# Fix CentOS 7 EOL (https://serverfault.com/a/1161847)
sed -i s/mirror.centos.org/vault.centos.org/g /etc/yum.repos.d/CentOS-*.repo
sed -i s/^#.*baseurl=http/baseurl=http/g /etc/yum.repos.d/CentOS-*.repo
sed -i s/^mirrorlist=http/#mirrorlist=http/g /etc/yum.repos.d/CentOS-*.repo

# Install C++17 Compiler
yum -y install devtoolset-9

# Build
PARAVIEW_DIR=\`find /builds/gitlab-kitware-sciviz-ci/build/install/lib/cmake -type d -regex "/builds/gitlab-kitware-sciviz-ci/build/install/lib/cmake/paraview-[0-9,.]*"\`
scl enable devtoolset-9 -- cmake -S /vofflow/src -B /vofflow/build -DParaView_DIR=\$PARAVIEW_DIR -DCMAKE_BUILD_TYPE=Release
scl enable devtoolset-9 -- cmake --build /vofflow/build --parallel 4
scl enable devtoolset-9 -- cmake --install /vofflow/build --prefix /plugin_build

# Cleanup static libs
find /plugin_build -name "*.a" -type f -delete

# Cleanup lib64 > lib
mv /plugin_build/lib64 /plugin_build/lib

# Cleanup file owners
chown -R $(id -u):$(id -g) /plugin_build
EOF
chmod +x "plugin_build/run.sh"

# Execute build in docker
docker run --rm -v "$(pwd)/plugin_build:/plugin_build" kitware/paraview_org-plugin-devel:5.13.2 /bin/bash -c "/plugin_build/run.sh"

# Cleanup
rm "plugin_build/run.sh"
