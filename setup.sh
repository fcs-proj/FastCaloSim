#!/bin/bash

# Use environment variable if set, otherwise use default version
KEY4HEP_VERSION=${KEY4HEP_VERSION:-"2024-11-28"}

# Setup key4hep environment
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${KEY4HEP_VERSION}"

# Set up python virtual environment only if it doesn't exist
mkdir -p ~/.local/bin
export PATH=~/.local/bin:"${PATH}"
curl -sL https://raw.githubusercontent.com/jbeirer/cvmfs-venv/main/cvmfs-venv.sh -o ~/.local/bin/cvmfs-venv
chmod +x ~/.local/bin/cvmfs-venv

# Check if the virtual environment already exists
if [ ! -d "fcs-venv" ]; then
    cvmfs-venv fcs-venv
else
    echo "Virtual environment 'fcs-venv' already exists. Skipping creation."
fi

# Activate the environment
. fcs-venv/bin/activate
