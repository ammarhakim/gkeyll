# The 2024 CUDA versions appear to be incompatible with GPU parallelism, loading the 2023 versions.
module load 2023
module load CUDA/12.1.1
module load OpenMPI/4.1.5-NVHPC-23.7-CUDA-12.1.1
module load NCCL/2.18.3-GCCcore-12.3.0-CUDA-12.1.1

cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --build-cudss=yes --prefix=$PREFIX