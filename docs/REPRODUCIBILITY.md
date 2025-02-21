# Reproducibility Steps for Visualization of Finite-Time Separation in Multiphase Flow

We tested this setup with Ubuntu 22.04.

## Content

- [Setup ParaView and the VofFlow Plugin](#setup-paraview-and-the-vofflow-plugin)
- [Download Dataset](#download-dataset)
- [Reproducing the Results](#reproducing-the-results)

## Setup ParaView and the VofFlow Plugin

### Option 1: Using Prebuild Binaries

Plugin binaries for ParaView 5.13.2 are provided on the [Release Page](https://github.com/UniStuttgart-VISUS/vof-flow/releases).
They are compatible with the official ParaView 5.13.2 MPI release for Linux available from the [ParaView website](https://www.paraview.org/download/):

- [Download ParaView 5.13.2 MPI Linux](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.13&type=binary&os=Linux&downloadFile=ParaView-5.13.2-MPI-Linux-Python3.10-x86_64.tar.gz)

Please unpack all downloaded files to `~/vofflow/install`, i.e., the `paraview` binary should be located at `~/vofflow/install/bin/paraview`.
This path is assumed in our example scripts.

> The released plugin binaries for Linux are compiled using the [GitHub Actions workflow in this repository](../.github/workflows/build.yml).
> The workflow is based on the [ParaViewEasyPluginBuilder](https://gitlab.kitware.com/paraview/paraview-easy-plugin-builder).
> The Windows version was compiled manually using the [ParaView Plugin Windows Binary Compatible Guide](https://gitlab.kitware.com/paraview/paraview-plugin-windows-binary-compatible-guide).

### Option 2: Build Plugin using the ParaViewEasyPluginBuilder

This option also uses the official ParaView release but compiles the plugin locally.
Please follow the steps from Option 1 to download and unpack the ParaView 5.13.2 release.
Please make sure Docker is installed on the system, see [https://docs.docker.com/engine/install/ubuntu/](https://docs.docker.com/engine/install/ubuntu/).

Then, run [`./build_plugin_linux.sh`](scripts/build_plugin_linux.sh) to build the plugin binaries.
The plugin will be created in a folder called `plugin_build` in the current directory.
Please move the content from `plugin_build` to `~/vofflow/install`.

> The script is based on the [ParaViewEasyPluginBuilder](https://gitlab.kitware.com/paraview/paraview-easy-plugin-builder).
> This is a Docker-based script to build ParaView plugins compatible with the official binary release of ParaView and does not require to build ParaView manually.
> We slightly modified the ParaViewEasyPluginBuilder by installing `devtoolset-9` for a C++17 compatible compiler.

### Option 3: Compile ParaView and the Plugin (Ubuntu 22.04)

- Install build environment:
  ```shell
  sudo apt install build-essential git git-lfs cmake ninja-build
  ```
- Install ParaView dependencies:
  ```shell
  sudo apt install qtbase5-dev qttools5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools libgl1-mesa-dev python3-dev libopenmpi-dev
  ```
- Install vcpkg dependencies:
  ```shell
  sudo apt install curl zip unzip tar pkg-config
  ```
- Clone ParaView and VofFlow:
  ```shell
  mkdir ~/vofflow && cd ~/vofflow
  git clone --recursive --branch v5.13.2 https://gitlab.kitware.com/paraview/paraview.git
  git clone https://github.com/UniStuttgart-VISUS/vof-flow.git
  ```
- Configure, Build, and Install ParaView:
  ```shell
  mkdir pv-build && cd pv-build
  cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DCMAKE_INSTALL_PREFIX=../install ../paraview
  cmake --build .
  cmake --install .
  cd ..
  ```
- Configure, Build, and Install VofFlow:
  ```shell
  mkdir vf-build && cd vf-build
  cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=~/vofflow/install -DCMAKE_INSTALL_PREFIX=../install ../vof-flow/
  cmake --build .
  cmake --install .
  cd ..
  ```

> This manual setup was tested with ParaView 5.13.2.

#### Alternative using system packages instead of vcpkg

By default, all dependencies except ParaView itself are downloaded automatically using vcpkg during the CMake configure of our project.
If you prefer system packages over vcpkg, you can install them manually:

```shell
sudo apt install libcgal-dev
```

In the VofFlow configure step, pass the additional parameter `-DVOFFLOW_USE_VCPKG=OFF`:

```shell
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=~/vofflow/install -DCMAKE_INSTALL_PREFIX=../install -DVOFFLOW_USE_VCPKG=OFF ../vof-flow/
```

#### Alternative for the ParaView-Superbuild

If the ParaView-Superbuild is used to build ParaView, there could be version conflicts between the Boost version of the ParaView-Superbuild and the vcpkg provided Boost version.
We provide a CMake option to disable downloading Boost with vcpkg: `-DVOFFLOW_DISABLE_VCPKG_BOOST=ON`.

## Download Dataset

We have published the jet-collision dataset from our paper here: [https://doi.org/10.18419/darus-4225](https://doi.org/10.18419/darus-4225).
We provide the script `download_data.py` to download all files automatically.
The script requires the Python package `requests`, e.g., run `pip install -r requirements.txt`.

The script can be run with different options to download the different versions of the dataset:

- `./download_data.py --ds 0` for the full-resolution dataset
- `./download_data.py --ds 1` for the downsampled variant
- `./download_data.py --ds 2` for the two times downsampled variant

Please move the downloaded datasets to `~/vofflow/data`.

## Reproducing the Results

### Prerequisites

Our reproducibility scripts assume that:

- The ParaView binary is installed at `~/vofflow/install/bin/paraview` (and `pvbatch`, etc.)
- The Plugin is installed at `~/vofflow/install/lib/paraview-5.13/plugins/VofFlow/VofFlow.so`
- The dataset "ds2" was downloaded at `~/vofflow/data/jet-collision-ds2/jet-collision-ds2.pvd`

All scripts have the paths configured at the top.
If your setup varies, you can update the paths accordingly, i.e., if you use one of the other datasets.

### Run Separation Boundaries

This step runs our algorithm and writes all resulting data, i.e., the separation boundaries, to disk in the VTK data formats.
The output is written to `~/vofflow/output/result`.

```shell
~/vofflow/install/bin/pvbatch separation_boundaries.py
```

Due to the computation time, we recommend running this step with MPI:

```shell
mpirun -np 32 ~/vofflow/install/bin/pvbatch separation_boundaries.py
```

### Run Surface Smoothing Post-Processing Step

The surface smoothing filter is not MPI aware.
It is recommended to run the post-processing step without MPI to avoid gaps at process boundaries.

```shell
~/vofflow/install/bin/pvbatch surface_smoothing.py
```

### Run Separation Boundary Rendering

Generates one timestep from Figure 14.

```shell
~/vofflow/install/bin/pvbatch figure_jet_bounds.py
```

### Run PLIC Rendering

Generates the PLIC surfaces of the jet shown in Figure 13.

```shell
~/vofflow/install/bin/pvbatch figure_jet_plic.py
```
