module unload darshan
: "${PREFIX:=/global/homes/m/mana/perlmutter/gkeyll/code/gkeyll_cpu/gkylsoft}"
./configure CC=cc --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --use-lua=yes;


