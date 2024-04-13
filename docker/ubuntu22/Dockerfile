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

RUN git clone --branch "v6-30-06" --depth=1 https://github.com/root-project/root.git root_src \
&& mkdir root_build \
&& cd root_build \
&& cmake -DCMAKE_INSTALL_PREFIX=../root ../root_src -D mathmore=ON \
&& cmake --build . --target install -j4 \
&& rm -rf ../root_src ../root_build

# Install lwtnn for inference
# See https://github.com/lwtnn/lwtnn for instructions
RUN git clone --branch "v2.13" --depth=1 https://github.com/lwtnn/lwtnn.git lwtnn_src \
 && mkdir lwtnn_build lwtnn_install\
 && cd lwtnn_build \
 && cmake -DBUILTIN_BOOST=true -DBUILTIN_EIGEN=true -DCMAKE_INSTALL_PREFIX=../lwtnn ../lwtnn_src \
 && make -j4 \
 && make install \
 && rm -rf ../lwtnn_src ../lwtnn_build

# Set ROOT environment
ENTRYPOINT ["/bin/bash", "-c", "source /root/bin/thisroot.sh;/bin/bash"]
