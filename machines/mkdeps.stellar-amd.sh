module load nvhpc/25.5
module load openmpi/nvhpc-25.5/4.1.8
module load openblas/0.3.x

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --prefix=$PREFIX --build-cudss=yes --build-luajit=yes
