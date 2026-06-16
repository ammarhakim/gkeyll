module load PrgEnv-gnu/8.6.0
module load craype-accel-nvidia80
module load cray-mpich/8.1.30
module load cudatoolkit/12.9
module load nccl/2.18.3-cu12
module load cray-libsci/25.09.0

: "${PREFIX:=/pscratch/sd/m/mana/gkeyll/code/gkeyll_gpu1/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --build-cudss=yes --prefix=$PREFIX --build-luajit=yes MPICC=mpicc  MPICXX=mpicxx
