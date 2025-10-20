# The 2024 CUDA versions appear to be incompatible with GPU parallelism, loading the 2023 versions.
echo "RECOMMENDED: add the following module loads to ~/.bashrc: 2023, CUDA/12.1.1, OpenMPI/4.1.5-NVHPC-23.7-CUDA-12.1.1, NCCL/2.18.3-GCCcore-12.3.0-CUDA-12.1.1"

module load 2023
module load CUDA/12.1.1
module load OpenMPI/4.1.5-NVHPC-23.7-CUDA-12.1.1
module load NCCL/2.18.3-GCCcore-12.3.0-CUDA-12.1.1
: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=nvcc ARCH_FLAGS="-march=native" CUDA_ARCH=80 --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --use-mpi=yes --mpi-inc=$EBROOTOPENMPI/include --mpi-lib=$EBROOTOPENMPI/lib64 --use-nccl=yes --nccl-inc=$EBROOTNCCL/include --nccl-lib=$EBROOTNCCL/lib64 --use-lua=yes --use-cudss=yes;