FROM ubuntu:22.04

ENV LANG=C.UTF-8

# Install packages listes in packages.txt
COPY packages.txt packages.txt

RUN apt-get update -qq \
 && ln -sf /usr/share/zoneinfo/UTC /etc/localtime  \
 && apt-get -y install $(cat packages.txt) --no-install-recommends \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && rm -rf /var/cache/apt/archives/* \
 && rm -rf /var/lib/apt/lists/* \
 && rm packages.txt

# See https://github.com/root-project/root-ci-images/blob/main/alma9/packages.txt for ROOT dependencies
# See https://root.cern/install/all_releases/ for available ROOT binaries
# We compile ROOT from source as currently we require mathmore which is off by default in the precompiled binaries
RUN git clone --branch "v6-30-06" --depth=1 https://github.com/root-project/root.git root_src \
&& mkdir root_build root_install \
&& cd root_build \
&& cmake -DCMAKE_INSTALL_PREFIX=../root-install ../root_src -D mathmore=ON \
&& cmake --build . --target install -j4 \
&& rm -rf ../root_src ../root_build

# Install lwtnn for inference
# See https://github.com/lwtnn/lwtnn for instructions
RUN git clone --branch "v2.13" --depth=1 https://github.com/lwtnn/lwtnn.git lwtnn_src \
 && mkdir lwtnn_build lwtnn_install\
 && cd lwtnn_build \
 && cmake -DBUILTIN_BOOST=true -DBUILTIN_EIGEN=true ../lwtnn_src \
 && make -j4 \
 && make install \
 && rm -rf ../lwtnn_src ../lwtnn_build

 # Install onnxruntime for inference (building onnx needs cmake-3.27 or higher)
 # See https://onnxruntime.ai/docs/build/inferencing.html for build instructions
 RUN git clone --branch "v1.16.3" --recursive https://github.com/Microsoft/onnxruntime.git onnx_src \
 && python3 -m ensurepip --default-pip && python3 -m pip install cmake \
 && mkdir onnx_build && cd onnx_build \
 && ../onnx_src/build.sh --allow_running_as_root --config RelWithDebInfo --build_shared_lib --parallel 4 --skip_submodule_sync --skip_tests --build_dir . \
 && make install -C RelWithDebInfo \
 && rm -rf ../onnx_src ../onnx_build


# Add ROOT .bashrc to the image 
# Note: this is just for the .devcontainer which ignores the docker
# jobs in CI will be launched in separate sub shells, ignoring any set environment variables
RUN touch /root/.bashrc \
 && echo 'source /root-install/bin/thisroot.sh' >> /root/.bashrc
