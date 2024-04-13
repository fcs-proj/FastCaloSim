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
 && mkdir root-build root-install \
 && cd root-build \
 && cmake -DCMAKE_INSTALL_PREFIX=../root-install ../root_src -D mathmore=ON \
 && cmake --build . --target install -j4


# Set ROOT environment
ENTRYPOINT ["/bin/bash", "-c", "source /root-install/bin/thisroot.sh;/bin/bash"]
