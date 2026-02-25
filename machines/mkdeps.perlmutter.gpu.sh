module load PrgEnv-gnu/8.6.0
module load craype-accel-nvidia80
module load cray-mpich/9.0.1
module load cudatoolkit/13.0
module load nccl/2.29.2-cu13
module load cray-libsci/25.09.0

: "${PREFIX:=/pscratch/sd/m/mana/gkeyll/code/gkeyll_gpu1/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --build-cudss=yes --prefix=$PREFIX --build-luajit=yes MPICC=mpicc  MPICXX=mpicxx --build-adas=yes
