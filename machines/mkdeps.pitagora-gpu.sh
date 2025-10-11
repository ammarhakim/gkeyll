# Loading modules, pitagora currently only supports cuda 12.6 

module load cuda/12.6
module load module load openmpi/4.1.6--gcc--12.3.0
module load nccl/2.22.3-1--gcc--12.3.0
module load cmake/3.27.9

cd install-deps
: "${PREFIX:=$HOME/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --build-cudss=yes --build-adas=yes --prefix=$PREFIX
