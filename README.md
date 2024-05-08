# Legendre-to-Chebyshev

[![L2C-CI](https://github.com/spectralDNS/Legendre-to-Chebyshev/actions/workflows/l2c.yml/badge.svg)](https://github.com/spectralDNS/Legendre-to-Chebyshev/actions/workflows/l2c.yml)

Routines for Legendre to Chebyshev (and inverse) transforms.

This is first and foremost the implementation of a Fast Multipole Method similar to (but not exactly like) the one described in

  * B. K. Alpert and V. Rokhlin, A fast algorithm for the evaluation of legendre expansions, 389 SIAM Journal on Scientific and Statistical Computing, 12 (1991), pp. 158–179, https://doi.390org/10.1137/0912009.391

The implemented method is described in the preprint [A faster Legendre-to-Chebyshev transform](https://github.com/spectralDNS/Legendre-to-Chebyshev/FMM_paper.pdf).

There are several implementations in the src directory:
  * python - A short, vectorized Python implementation
  * C - An efficient C implementation
  * cython - Wrapping of the fast C code using cython
  * ctypes - Wrapping of the fast C code using ctypes
  * Multiprec - Multiprecision implementations of a direct method using both Python and C++. The multiprecision methods are only used for verification.

In addition there is
  * bin - An executable l2c that can run various tests

# Installation
The code is set up to be compiled with the [meson](https://mesonbuild.com) build system. It should work by cloning this repository and then

    cd src
    meson setup build
    meson install -C build

If you want to install just locally, then use, e.g.,

    meson setup build --prefix=$PWD/build-install

The installation can be tested with

    meson test -C build

# Codespace
Another simple way to test this code is to create a codespace. The environment.yml file in the root folder will then make sure that the codespace creates a conda environment with all necessary dependencies already built. Just press the codespace button and wait awhile for the environment to build. Then enable the environment and run some tests or test the executable `l2c`

     source activate ./venv
     cd src
     meson setup build --prefix=$PWD/build-install --includedir=$CONDA_PREFIX/include --libdir=$CONDA_PREFIX/lib
     meson install -C build
     export PATH=$PWD/build-install/bin:$PATH
     l2c -N1000 -d2 # runs a forward and backward transform and computes the error for an array of length 1000
     meson test -C build # runs all the tests
