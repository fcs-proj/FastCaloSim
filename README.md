# FastCaloSim
 [![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](code_of_conduct.md) 

FastCaloSim is an experiment-independent toolkit for the fast parametrised and ML-based simulation of electromagnetic and hadronic showers in (high energy) physics experiments implemented in C++.

More information can be found in the [documentation](https://fastcalosim.readthedocs.io/) and the [doxygen](https://jbeirer.github.io/FastCaloSim/index.html).

# Quick Start

FastCaloSim is developed in C++ and is build using [CMake](https://cmake.org). The
following commands will clone the repository, configure, and build the library

```sh
git clone https://github.com/fcs-project/FastCaloSim <source>
cmake -B <build> -S <source> -D CMAKE_BUILD_TYPE=Release
cmake --build <build>
```
For install options and instruction on how to include FastCaloSim in your experiment see the [BUILDING](BUILDING.md) document. For advanced developer configuration with [presets][1] and other useful information see the [HACKING](HACKING.md) document.

# Contributing

See the [CONTRIBUTING](CONTRIBUTING.md) document.

# Authors

Contributors to the FastCaloSim project are listed in the [AUTHORS](AUTHORS) file.

# Licensing

The FastCaloSim project is published under the terms of the Apache 2.0 license, except where other licenses apply. A copy of the license can be found in the in the [LICENSE](LICENSE) file or at https://www.apache.org/licenses/LICENSE-2.0



[1]: https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html
