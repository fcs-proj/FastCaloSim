FROM gitlab-registry.cern.ch/linuxsupport/alma9-base

# Install packages listes in packages.txt
COPY packages.txt packages.txt

RUN dnf update -y \
 && dnf install -y epel-release \
 && dnf install -y 'dnf-command(config-manager)' \
 && crb enable \
 && dnf install -y --setopt=install_weak_deps=False $(cat packages.txt) \
 && rm -f packages.txt \
 && dnf clean all \
 && rm -rf /var/cache/dnf

# Upgrade CMake version
# TODO: delete .sh file after installation
ARG CMAKE_VERSION=3.29.2
RUN wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.sh \
 && chmod +x cmake-${CMAKE_VERSION}-linux-x86_64.sh \  
 && sh cmake-${CMAKE_VERSION}-linux-x86_64.sh --skip-license --prefix=/usr/local

# Build ROOT from source.
# Note: mathmore is currently needed for interpolation, clad is deactivated due to build issues
# Important to install it in /usr/local to make it available in the PATH
ARG ROOT_VERSION=6-30-06
RUN git clone --branch v${ROOT_VERSION} --depth=1 https://github.com/root-project/root.git root_src \
 && mkdir root_build \
 && cd root_build \
 && cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../root_src -D mathmore=ON -D clad=OFF \
 && cmake --build . --target install -j6 \
 && rm -rf ../root_src ../root_build


# Build onnxruntime for inference from source
# See https://onnxruntime.ai/docs/build/inferencing.html for build instructions
ARG ONNXRUNTIME_VERSION=1.16.3
RUN git clone --branch v${ONNXRUNTIME_VERSION} --recursive https://github.com/Microsoft/onnxruntime.git onnx_src \
 && mkdir onnx_build && cd onnx_build \
 && ../onnx_src/build.sh --allow_running_as_root --config RelWithDebInfo --build_shared_lib --parallel 4 --skip_submodule_sync --skip_tests --build_dir . \
 && make install -C RelWithDebInfo \
 && rm -rf ../onnx_src ../onnx_build

# Build lwtnn for inference from source
# See https://github.com/lwtnn/lwtnn for instructions
ARG LWTNN_VERSION=2.13
RUN git clone --branch v${LWTNN_VERSION} --depth=1 https://github.com/lwtnn/lwtnn.git lwtnn_src \
 && mkdir lwtnn_build lwtnn_install\
 && cd lwtnn_build \
 && cmake -DBUILTIN_BOOST=true -DBUILTIN_EIGEN=true -DCMAKE_INSTALL_PREFIX=/usr/local ../lwtnn_src \
 && make -j4 \
 && make install \
 && rm -rf ../lwtnn_src ../lwtnn_build