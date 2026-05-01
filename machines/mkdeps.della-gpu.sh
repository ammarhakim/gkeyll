module load nvhpc/25.5 # For NCCL.
module load cudatoolkit/12.9
module load openmpi/cuda-12.9/nvhpc-25.5/4.1.8

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes --prefix=$PREFIX --build-cudss=yes --build-luajit=yes
