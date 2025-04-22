# FastCaloSim
 [![Checks](https://github.com/fcs-proj/FastCaloSim/actions/workflows/pipeline.yml/badge.svg)](https://github.com/fcs-proj/FastCaloSim/actions/workflows/pipeline.yml)
 [![Release](https://img.shields.io/github/v/release/fcs-proj/fastcalosim?include_prereleases)](https://github.com/fcs-proj/fastcalosim/releases)
 [![Coverage](https://codecov.io/gh/fcs-proj/FastCaloSim/graph/badge.svg)](https://codecov.io/gh/fcs-proj/FastCaloSim)
 [![Doxygen Documentation](https://img.shields.io/badge/docs-Doxygen-blue)](https://fcs-proj.github.io/FastCaloSim/)
 [![DOI](https://zenodo.org/badge/781604786.svg)](https://doi.org/10.5281/zenodo.14308464)

FastCaloSim is an experiment-independent toolkit for the fast parametrised and ML-based simulation of electromagnetic and hadronic showers in (high energy) physics experiments implemented in C++.

# Quick Start

FastCaloSim is developed in C++ and is build using [CMake](https://cmake.org). The
following commands will clone the repository, configure, and build the library

```sh
git clone https://github.com/fcs-proj/FastCaloSim <source>
cmake -B <build> -S <source> -D CMAKE_BUILD_TYPE=Release
cmake --build <build>
```
For install options and instruction on how to include FastCaloSim in your experiment see the [BUILDING](BUILDING.md) document. For advanced developer configuration with [presets][1] and other useful information see the [HACKING](HACKING.md) document.

# Development in Containers

We recommend using Docker containers for FastCaloSim development. We provide two official containers with all necessary dependencies pre-installed:

- AlmaLinux 9 based: [`fcsproj/fastcalosim-alma9`](https://hub.docker.com/r/fcsproj/fastcalosim-alma9/tags)
- Ubuntu 24.04 based: [`fcsproj/ubuntu24-fastcalosim`](https://hub.docker.com/r/fcsproj/ubuntu24-fastcalosim/tags)

Both containers support x86_64 (Intel/AMD) and ARM64 (Apple Silicon) architectures, ensuring consistent development across platforms.

## VS Code Integration

When using VS Code with the DevContainer extension, containers are automatically pulled and configured when you open the repository, providing a complete development environment with all tools pre-configured.

## LXPLUS Development with Apptainer

For development on LXPLUS, you can use our Docker containers via Apptainer. We provide an `apptainer.sh` script in the repository that contains the necessary command:

```bash
# Run using the provided script
source apptainer.sh

# Or use the command directly
singularity exec --contain --cleanenv \
  -B $(pwd):/FastCaloSim \
  -B /tmp \
  docker://fcsproj/fastcalosim-alma9:main \
  bash
```
 
# Contributing

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

See the [CONTRIBUTING](CONTRIBUTING.md) document.

# Authors

Contributors to the FastCaloSim project are listed in the [AUTHORS](AUTHORS) file.

# Licensing

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

The FastCaloSim project is published under the terms of the Apache 2.0 license, except where other licenses apply. A copy of the license can be found in the in the [LICENSE](LICENSE) file or at https://www.apache.org/licenses/LICENSE-2.0



[1]: https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html
