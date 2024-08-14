#!/bin/bash

#export LD_LIBRARY_PATH=$PWD/build-install/lib  
#export DYLD_LIBRARY_PATH=$PWD/build-install/lib  
#export PYTHONPATH=$PYTHONPATH:$PWD/build-install/lib/python3.12/site-packages

export USE_ACCELERATE=0

rm -rf build
rm -rf build-install

if [ "$USE_ACCELERATE" -eq 1 ]; then
  export PATH="/opt/homebrew/opt/cython/bin:$PATH"
  meson setup build --prefix=$PWD/build-install --includedir=/opt/homebrew/include --libdir=/opt/homebrew/lib
else
  meson setup build --prefix=$PWD/build-install --includedir=$CONDA_PREFIX/include --libdir=$CONDA_PREFIX/lib
fi
meson configure -Dbuildtype=release -Doptimization=3 -Dc_args="-march=native -Ofast -fPIC" -Dc_link_args="-lm -fPIC"  build
meson compile -v -C build
meson install -C build
