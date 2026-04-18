# Copyright (c) 2026 CERN for the benefit of the FastCaloSim project
FROM gitlab-registry.cern.ch/linuxsupport/alma10-base

# Build type for all dependencies
ARG BUILD_TYPE=RelWithDebInfo

# Install extra packages and ROOT dependencies
# Note: ROOT packages taken from https://github.com/root-project/root-ci-images/blob/main/alma10/packages.txt
COPY extra_packages.txt extra_packages.txt
COPY root_dependencies.txt root_dependencies.txt

RUN dnf update -y \
 && dnf install -y epel-release \
 && dnf install -y 'dnf-command(config-manager)' \
 && crb enable \
 && dnf install -y --setopt=install_weak_deps=False \
    $(cat extra_packages.txt) $(cat root_dependencies.txt) \
 && dnf clean all \
 && rm -rf extra_packages.txt root_dependencies.txt /var/cache/dnf

# Build CMake from source
# See https://cmake.org/install/ for build instructions
ARG CMAKE_VERSION=3.29.2
RUN wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}.tar.gz \
 && tar -xzvf cmake-${CMAKE_VERSION}.tar.gz \
 && cd cmake-${CMAKE_VERSION} \
 && ./bootstrap --parallel=$(nproc) --prefix=/usr/local \
 && make -j$(nproc) \
 && make install \
 && cd .. \
 && rm -rf cmake-${CMAKE_VERSION} cmake-${CMAKE_VERSION}.tar.gz

# Build ROOT from source
# Note: mathmore is needed for interpolation, clad is deactivated due to build issues
# See https://root.cern/install/build_from_source/ for build instructions
ARG ROOT_VERSION=6-30-06
RUN git clone --branch v${ROOT_VERSION} --depth=1 https://github.com/root-project/root.git root_src \
 && mkdir root_build && cd root_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -Dmathmore=ON \
          -Dclad=OFF \
          ../root_src \
 && cmake --build . --target install -j$(nproc) \
 && rm -rf ../root_build

# Build Geant4 from source with GDML support and installed data
# Note: -ffp-contract=off -fno-fast-math ensures floating point reproducibility across optimization levels
# See https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/buildtools.html for build instructions
ARG GEANT4_VERSION=11.2.1
RUN git clone --depth=1 --branch=v${GEANT4_VERSION} https://github.com/Geant4/geant4 geant4_src \
 && mkdir geant4_build && cd geant4_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DCMAKE_CXX_FLAGS="-ffp-contract=off -fno-fast-math" \
          -DGEANT4_USE_GDML=ON \
          -DGEANT4_INSTALL_DATA=ON \
          -DGEANT4_INSTALL_DATA_TIMEOUT=3600 \
          ../geant4_src \
 && cmake --build . --target install -j$(nproc) \
 && rm -rf ../geant4_build

# Build onnxruntime for inference from source
# See https://onnxruntime.ai/docs/build/inferencing.html for build instructions
ARG ONNXRUNTIME_VERSION=1.24.4
RUN git clone --branch v${ONNXRUNTIME_VERSION} --recursive https://github.com/Microsoft/onnxruntime.git onnx_src \
 && mkdir onnx_build && cd onnx_build \
 && ../onnx_src/build.sh \
    --allow_running_as_root \
    --config ${BUILD_TYPE} \
    --build_shared_lib \
    --parallel $(nproc) \
    --compile_no_warning_as_error \
    --skip_submodule_sync \
    --skip_tests \
    --build_dir . \
    --cmake_extra_defines CMAKE_INSTALL_PREFIX=/usr/local \
 && cd ${BUILD_TYPE} \
 && cmake --build . --target install -j$(nproc) \
 && rm -rf ../../onnx_build

# Build lwtnn for inference from source
# See https://github.com/lwtnn/lwtnn for instructions
ARG LWTNN_VERSION=2.13
RUN git clone --branch v${LWTNN_VERSION} --depth=1 https://github.com/lwtnn/lwtnn.git lwtnn_src \
 && mkdir lwtnn_build && cd lwtnn_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DBUILTIN_BOOST=ON \
          -DBUILTIN_EIGEN=ON \
          ../lwtnn_src \
 && cmake --build . --target install -j$(nproc) \
 && rm -rf ../lwtnn_build

# Build CLHEP from source
# See https://gitlab.cern.ch/CLHEP/CLHEP for instructions
ARG CLHEP_VERSION=2_4_7_1
RUN git clone --branch CLHEP_${CLHEP_VERSION} --depth=1 https://gitlab.cern.ch/CLHEP/CLHEP.git clhep_src \
 && mkdir clhep_build && cd clhep_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          ../clhep_src \
 && cmake --build . --target install -j$(nproc) \
 && rm -rf ../clhep_build

# Build HDF5 from source with C++ and high-level bindings
# See https://github.com/HDFGroup/hdf5 for instructions
ARG HDF5_VERSION=1.14.6
RUN git clone --branch hdf5_${HDF5_VERSION} --depth=1 https://github.com/HDFGroup/hdf5.git hdf5_src \
 && mkdir hdf5_build && cd hdf5_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local \
          -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
          -DBUILD_STATIC_LIBS=OFF \
          -DHDF5_BUILD_CPP_LIB=ON \
          -DHDF5_BUILD_TOOLS=OFF \
          -DHDF5_BUILD_EXAMPLES=OFF \
          -DBUILD_TESTING=OFF \
          ../hdf5_src \
 && cmake --build . --target install -j$(nproc) \
 && rm -rf ../hdf5_build

# Create a virtual environment
ENV VIRT_ENV=/opt/venv
RUN python3 -m venv $VIRT_ENV
ENV PATH="$VIRT_ENV/bin:$PATH"

# Install python packages
# Note: pygeosimplify installed without dependencies due to https://github.com/g4edge/pyg4ometry/issues/243
COPY python_packages.txt python_packages.txt
RUN pip3 install --upgrade pip \
 && pip3 install $(cat python_packages.txt) \
 && pip3 install pygeosimplify --no-deps \
 && rm -rf python_packages.txt
