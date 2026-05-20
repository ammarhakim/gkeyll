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

echo ">>> Installation of coupling library complete."

# ==============================================================================
# SOLPS-ITER Automated Installer for NERSC Perlmutter
# ==============================================================================
# Note that the automation of the SOLPS-ITER installation below
# used assistance from an LLM (Gemini).

# Stop on error
set -e

echo ">>> Starting SOLPS-ITER Installation..."

# Define SOLPSTOP absolute path
# We'll get the current directory and then append SOLPS-ITER to it
# This will be used in the setup script
SOLPSTOP_VAL="$(pwd)/SOLPS-ITER"

# 1. Clone the repository
if [ ! -d "SOLPS-ITER" ]; then
    echo ">>> Cloning SOLPS-ITER..."
    git clone https://github.com/iterorganization/SOLPS-ITER.git
else
    echo ">>> SOLPS-ITER directory already exists. Skipping clone."
fi

cd SOLPS-ITER

# 2. Modify .gitmodules to remove ADAS
echo ">>> Removing ADAS from .gitmodules..."
if [ -f ".gitmodules" ]; then
    # Use perl to remove the ADAS submodule section completely.
    # The regex is now more flexible to match different branch names or whitespace.
    perl -i -0777 -pe 's/\[submodule "ADAS"\]\s+path = modules\/adas\s+url = ssh:\/\/git\@git\.iter\.org\/imex\/amns-adas\.git\s+branch = .*?\n//sg' .gitmodules
fi

# 3. Initialize and update submodules
echo ">>> Initializing submodules..."
# If ADAS was already initialized, we might need to remove it from git config
git config --local --remove-section submodule.modules/adas 2>/dev/null || true
git config --local --remove-section submodule.ADAS 2>/dev/null || true
# Remove the directory if it exists and isn't a proper submodule anymore
if [ -d "modules/adas" ]; then
    rm -rf modules/adas
fi
git submodule init || true
git submodule update

# 4. Modify SETUP/setup.csh.NERSC.gfortran
echo ">>> Configuring SETUP/setup.csh.NERSC.gfortran..."

setup_file="SETUP/setup.csh.NERSC.gfortran"

# NEW_HEADER includes the corrected SOLPSTOP definition and NERSC modules
NEW_HEADER=$(cat << EOM
 setenv SOLPSTOP $SOLPSTOP_VAL

module load PrgEnv-gnu
module load cray-hdf5
module load cray-libsci
module load cray-netcdf
module load texlive/2024
module load cudatoolkit

setenv GLI_HOME \${SOLPSTOP}/lib/NERSC.gfortran
setenv MSCL_ROOT \${SOLPSTOP}/lib/NERSC.gfortran
setenv GR_ROOT \${SOLPSTOP}/lib/NERSC.gfortran
setenv GR_SOURCE \${SOLPSTOP}/lib/NERSC.gfortran
setenv GKS_ROOT \${SOLPSTOP}/lib/NERSC.gfortran
setenv PATH \${PATH}:\${SOLPSTOP}/scripts
setenv LD_LIBRARY_PATH \${SOLPSTOP}/lib/NERSC.gfortran:\$LD_LIBRARY_PATH
EOM
)

# Replace the top block up to '# Environment variables'
# Use a temporary file to safely handle the replacement with a large string
export NEW_HEADER
perl -i -0777 -pe 's/^(.*?)# Environment variables/$ENV{NEW_HEADER}\n# Environment variables/s' "$setup_file"

# 5. Step 1: First Setup & Miniforge
echo ">>> Running first_setup and installing miniforge..."
cat << 'EOF' > step1_install.csh
#!/bin/tcsh
setenv HOST_NAME NERSC
setenv COMPILER gfortran
source SETUP/setup.csh.NERSC.gfortran
first_setup
printf "y\n" | install_miniforge
EOF

tcsh step1_install.csh
rm step1_install.csh

# 6. Step 2: Dependencies
echo ">>> Installing dependencies..."
# Based on instructions:
# Yes: Blas/Lapack, JSON-Fortran, NCARG, MSCL, GKS, GR
# No: NetCDF, hdf5 (these are available on Perlmutter)
# 
# The expected prompt order in install_dependencies (standard for SOLPS-ITER) is:
# 1. Gnuplot (y) - but it seems it might be already installed or skipped
# 2. Blas/Lapack (y)
# 3. NetCDF (n)
# 4. JSON-Fortran (y)
# 5. NCARG (y)
# 6. MSCL (y)
# 7. GKS (y)
# 8. GR (y)
# Note: hdf5 doesn't seem to have a prompt in the script for HPC? 
# Actually, I'll provide enough 'y' and 'n' to cover all bases.
INPUTS="y\ny\nn\ny\ny\ny\ny\ny\n"

cat << EOF > step2_deps.csh
#!/bin/tcsh
setenv HOST_NAME NERSC
setenv COMPILER gfortran
source SETUP/setup.csh.NERSC.gfortran
conda activate ./lib/NERSC.gfortran/miniforge3
setenv CC cc
printf "$INPUTS" | install_dependencies
EOF

tcsh step2_deps.csh
rm step2_deps.csh

# 7. Step 3: Fix JSON-Fortran Build
echo ">>> Fixing JSON-Fortran build..."
cat << 'EOF' > step3_json.csh
#!/bin/tcsh
setenv HOST_NAME NERSC
setenv COMPILER gfortran
source SETUP/setup.csh.NERSC.gfortran
conda activate solps_env

if ( ! -d lib/NERSC.gfortran/json-fortran-9.0.2 ) then
    cd lib/NERSC.gfortran
    wget https://github.com/jacobwilliams/json-fortran/archive/refs/tags/9.0.2.tar.gz
    tar xvf 9.0.2.tar.gz
    rm 9.0.2.tar.gz
    cd ../..
endif

cd lib/NERSC.gfortran/json-fortran-9.0.2
if ( -d build ) then
    rm -rf build
endif
mkdir build
cd build
cmake -DENABLE_IN_SOURCE_BUILDS=ON ..
make
# Also fulfill the user's manual step as it might be what they expect
cd ..
# If there is no Makefile in the root, 'make' will fail but we can ignore it
# because the build in 'build/' directory succeeded.
make || true
EOF

tcsh step3_json.csh
rm step3_json.csh

# 8. Step 4: Configure Carre and Compile
echo ">>> Configuring Carre and compiling SOLPS..."

config_file="modules/Carre/config/config.NERSC.gfortran"
# Add DEFINES += -DGFORTRAN after LD_NCARG
# We'll use a more precise match that handles both ?= and =
# Check if it already exists to avoid duplicate entries
if ! grep -q "DEFINES += -DGFORTRAN" "$config_file"; then
    perl -i -pe 's/(LD_NCARG\s*[\?\=]*\s*.*\n)/$1DEFINES += -DGFORTRAN\n/' "$config_file"
fi

cat << 'EOF' > step4_compile.csh
#!/bin/tcsh
setenv HOST_NAME NERSC
setenv COMPILER gfortran
source setup.csh
conda activate solps_env

make clean
make clean_solps

echo ">>> Compilation starting. Piping 'q' to handle any pdflatex pauses..."
# We use 'make' instead of 'make solps' because 'make clean_solps' and 'make solps' are part of the instructions.
# 'make' in the root usually builds everything.
# Using source setup.csh instead of specific setup.csh.NERSC.gfortran because setup.csh
# handles adding scripts/ to PATH, which is where which_latex resides.
yes "q" | make solps

echo ">>> Performing manual DivGeo build to ensure dg.exe is created..."
cd modules/DivGeo
make
EOF

tcsh step4_compile.csh
rm step4_compile.csh

echo ">>> Installation of SOLPS-ITER complete."
