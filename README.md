# Legendre-to-Chebyshev

[![L2C-CI](https://github.com/spectralDNS/Legendre-to-Chebyshev/actions/workflows/l2c.yml/badge.svg)](https://github.com/spectralDNS/Legendre-to-Chebyshev/actions/workflows/l2c.yml)

Routines for Legendre to Chebyshev (and inverse) transforms.

This is first and foremost the implementation of a Fast Multipole Method similar to the one described in

  * B. K. Alpert and V. Rokhlin, A fast algorithm for the evaluation of legendre expansions, 389 SIAM Journal on Scientific and Statistical Computing, 12 (1991), pp. 158–179, https://doi.390 org/10.1137/0912009.391

There are several implementations in the src directory:
  * python - A short, vectorized Python implementation
  * C - An efficient C implementation
  * cython - Wrapping of the fast C code
  * Multiprec - Multiprecision implementations of a direct method using both Python and C++. The multiprecision methods are only used for verification.

In addition there are
  * test - Some test routines
  * bin - An executable that can run various tests

# Installation
The code is set up to be compiled with the [meson](https://mesonbuild.com) build system. It should work by cloning this repository and then 

    cd src
    meson setup build 
    meson install -C build
  
If you want to install just locally, then use, e.g.,

    meson setup build --prefix=$PWD/build-install

The installation can be tested with

    meson test -C build
    
