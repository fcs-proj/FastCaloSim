#!/usr/bin/env bash

# Script to set up python virtual environment on top of lcg or key4hep environment
set -e

mkdir -p ~/.local/bin
export PATH=~/.local/bin:"${PATH}"
curl -sL https://raw.githubusercontent.com/jbeirer/cvmfs-venv/main/cvmfs-venv.sh \
    -o ~/.local/bin/cvmfs-venv
chmod +x ~/.local/bin/cvmfs-venv
cvmfs-venv fcs-venv  && . fcs-venv/bin/activate
