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

# ROOT install from precompiled binaries
# See https://github.com/root-project/root-ci-images/blob/main/ubuntu22/packages.txt for ROOT dependencies
# See https://root.cern/install/all_releases/ for available ROOT binaries
ARG ROOT_BINARIES=root_v6.30.06.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz
RUN wget https://root.cern/download/$ROOT_BINARIES \
    && tar -xzvf $ROOT_BINARIES \
    && rm $ROOT_BINARIES

# Set ROOT environment
ENTRYPOINT ["/bin/bash", "-c", "source /root/bin/thisroot.sh;/bin/bash"]
