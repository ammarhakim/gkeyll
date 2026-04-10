#!/bin/bash

# ==============================================================================
# Gk-Neutral_coupling Automated Installer
# ==============================================================================

echo ">>> Starting GK-Neutral_coupling Installation..."
cd ../../gkylsoft/

if [ ! -d "GK-Neutral_coupling" ]; then
    echo ">>> Cloning GK-Neutral_coupling..."
    git clone https://github.com/jRoeltgen/GK-Neutral_coupling.git
else
    echo ">>> GK-Neutral_coupling directory already exists. Skipping clone."
fi

echo ">>> Installation Complete."
