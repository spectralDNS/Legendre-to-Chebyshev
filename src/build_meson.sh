#!/bin/bash

rm -rf build
rm -rf build-install
meson setup build --prefix=$PWD/build-install --includedir=$CONDA_PREFIX/include --libdir=$CONDA_PREFIX/lib
meson configure -Dbuildtype=release -Doptimization=3 -Dc_args="-march=native -Ofast -fPIC -fstack-protector-strong -pipe" -Dc_link_args="-lm -fPIC -fstack-protector-strong -pipe -D_FORTIFY_SOURCE=2"  build
meson compile -v -C build
meson install -C build