# Copyright (c) 2024 CERN for the benefit of the FastCaloSim project

FROM gitlab-registry.cern.ch/linuxsupport/alma9-base

ENV PIP_ROOT_USER_ACTION=ignore

# Install ROOT dependencies
# Note: Taken from https://github.com/root-project/root-ci-images/blob/main/alma9/packages.txt
COPY root_dependencies.txt root_dependencies.txt

RUN dnf update -y \
 && dnf install -y epel-release \
 && dnf install -y 'dnf-command(config-manager)' \
 && crb enable \
 && dnf install -y --skip-broken --setopt=install_weak_deps=False $(cat root_dependencies.txt) \
 && dnf clean all \
 && rm -rf root_dependencies.txt /var/cache/dnf

# Upgrade CMake version
ARG CMAKE_VERSION=3.29.2
RUN dnf install -y wget \
 && wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh \
 && chmod +x cmake-${CMAKE_VERSION}-linux-x86_64.sh \
 && sh cmake-${CMAKE_VERSION}-linux-x86_64.sh --skip-license --prefix=/usr/local \
 && rm cmake-${CMAKE_VERSION}-linux-x86_64.sh

# Build ROOT from source.
# Note: mathmore is currently needed for interpolation, clad is deactivated due to build issues
# Important to install it in /usr/local to make it available in the PATH
# See https://root.cern/install/build_from_source/ for build instructions
ARG ROOT_VERSION=6-30-06
RUN git clone --branch v${ROOT_VERSION} --depth=1 https://github.com/root-project/root.git root_src \
 && mkdir root_build \
 && cd root_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../root_src -D mathmore=ON -D clad=OFF \
 && cmake --build . --config RelWithDebInfo --target install -j6 \
 && rm -rf ../root_build

# Build Geant4 from source with GDML support and installed data
# See https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/buildtools.html for build instructions
ARG GEANT4_VERSION=11.2.1
RUN dnf install -y expat-devel xerces-c-devel \
 && git clone --depth=1 --branch=v${GEANT4_VERSION} https://github.com/Geant4/geant4 geant4_src \
 && mkdir geant4_build && cd geant4_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local -DGEANT4_USE_GDML=ON -DGEANT4_INSTALL_DATA=ON -DCMAKE_BUILD_TYPE=Debug ../geant4_src \
 && cmake --build . --config RelWithDebInfo --target install -j6  \
 && rm -rf ../geant4_build

# Build onnxruntime for inference from source
# See https://onnxruntime.ai/docs/build/inferencing.html for build instructions
ARG ONNXRUNTIME_VERSION=1.16.3
RUN git clone --branch v${ONNXRUNTIME_VERSION} --recursive https://github.com/Microsoft/onnxruntime.git onnx_src \
 && mkdir onnx_build && cd onnx_build \
 && ../onnx_src/build.sh --allow_running_as_root --config RelWithDebInfo --build_shared_lib --parallel 4 --skip_submodule_sync --skip_tests --build_dir . \
 && make install -C RelWithDebInfo \
 && rm -rf ../onnx_build

# Build lwtnn for inference from source
# See https://github.com/lwtnn/lwtnn for instructions
ARG LWTNN_VERSION=2.13
RUN git clone --branch v${LWTNN_VERSION} --depth=1 https://github.com/lwtnn/lwtnn.git lwtnn_src \
 && mkdir lwtnn_build && cd lwtnn_build \
 && cmake -DBUILTIN_BOOST=true -DBUILTIN_EIGEN=true -DCMAKE_INSTALL_PREFIX=/usr/local ../lwtnn_src \
 && make -j4 \
 && make install \
 && rm -rf ../lwtnn_build

# Build CLHEP from source
# See https://gitlab.cern.ch/CLHEP/CLHEP for instructions
ARG CLHEP_VERSION=2_4_7_1
RUN git clone https://gitlab.cern.ch/CLHEP/CLHEP.git --branch CLHEP_${CLHEP_VERSION} clhep_src \
 && mkdir clhep_build && cd clhep_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../clhep_src \
 && cmake --build . --config RelWithDebInfo --target install -j6 \
 && rm -rf ../clhep_build

# Install extra packages
COPY extra_packages.txt extra_packages.txt

RUN  dnf clean all && dnf update -y \
 && dnf install -y --setopt=install_weak_deps=False $(cat extra_packages.txt) \
 && dnf clean all \
 && rm -rf extra_packages.txt /var/cache/dnf

 # Install python packages
 COPY python_packages.txt python_packages.txt

 RUN python3 -m ensurepip --default-pip \
 && python3 -m pip install --upgrade pip \
 && pip3 install $(cat python_packages.txt) \
 && rm -rf python_packages.txt
