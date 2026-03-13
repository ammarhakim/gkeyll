module load nvhpc/25.5 # For NCCL.
module load cudatoolkit/12.9
module load openmpi/cuda-12.9/nvhpc-25.5/4.1.8
module load openblas/0.3.x

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --prefix=$PREFIX --build-cudss=yes --build-luajit=yes
