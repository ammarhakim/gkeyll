#!/bin/bash

# ==============================================================================
# Gk-Neutral_coupling Automated Installer
# ==============================================================================

echo ">>> Starting GK-Neutral_coupling Installation..."
cd ../../gkylsoft/

if [ ! -d "GK-Neutral_coupling" ]; then
    echo ">>> Downloading GK-Neutral_coupling..."
    wget https://github.com/jRoeltgen/GK-Neutral_coupling/archive/refs/tags/gkeyll.tar.gz
    mkdir GK-Neutral_coupling && tar -xzf gkeyll.tar.gz -C GK-Neutral_coupling --strip-components=1
    rm gkeyll.tar.gz
else
    echo ">>> GK-Neutral_coupling directory already exists. Skipping clone."
fi

echo ">>> Installation Complete."
