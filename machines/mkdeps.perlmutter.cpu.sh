module unload darshan
cd install-deps
: "${PREFIX:=/global/homes/m/mana/perlmutter/gkeyll/code/gkeyll_cpu/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --prefix=$PREFIX --build-adas=yes
