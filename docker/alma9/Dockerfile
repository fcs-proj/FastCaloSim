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


# ROOT install from precompiled binaries
# See https://github.com/root-project/root-ci-images/blob/main/alma9/packages.txt for ROOT dependencies
# See https://root.cern/install/all_releases/ for available ROOT binaries
ARG ROOT_BINARIES=root_v6.30.06.Linux-almalinux9.3-x86_64-gcc11.4.tar.gz
RUN wget https://root.cern/download/$ROOT_BINARIES \
    && tar -xzvf $ROOT_BINARIES \
    && rm $ROOT_BINARIES

# Set ROOT environment
ENTRYPOINT ["/bin/bash", "-c", "source /root/bin/thisroot.sh;/bin/bash"]