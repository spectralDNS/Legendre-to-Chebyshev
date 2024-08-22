# Legendre-to-Chebyshev

[![L2C-CI](https://github.com/spectralDNS/Legendre-to-Chebyshev/actions/workflows/l2c.yml/badge.svg)](https://github.com/spectralDNS/Legendre-to-Chebyshev/actions/workflows/l2c.yml)

Routines for Legendre to Chebyshev (and inverse) transforms.

This is first and foremost the implementation of a Fast Multipole Method similar to (but not exactly like) the one described in

  * B. K. Alpert and V. Rokhlin, A fast algorithm for the evaluation of legendre expansions, 389 SIAM Journal on Scientific and Statistical Computing, 12 (1991), pp. 158–179, https://doi.390org/10.1137/0912009.391

The implemented method is described in the preprint [A faster Legendre-to-Chebyshev transform](https://github.com/spectralDNS/Legendre-to-Chebyshev/blob/main/FMM_paper.pdf).

There are several implementations in the src directory:
  * python - A short, vectorized Python implementation
  * C - An efficient C implementation
  * cython - Wrapping of the fast C code using cython
  * ctypes - Wrapping of the fast C code using ctypes
  * Multiprec - Multiprecision implementations of a direct method using both Python and C++. The multiprecision methods are only used for verification.

In addition there is
  * bin - An executable l2c that can run various tests
  * results - Bash and Python files that can be used to recreate the figures in the paper above. They have obvious names, like table1.py. Note that the bash-files have precomputed data included and if you want to recompute something, then you need to open these files and modify at the top. For example, modify `rerun_fmm="no"` to `rerun_fmm="yes"` in order to rerun the results for the faster multipole method described in the paper.

# Installation
The code is set up to be compiled with the [meson](https://mesonbuild.com) build system. If all dependencies are easily found, then it should work by cloning this repository and then

    cd src
    meson setup build
    meson install -C build

If you want to install just locally, then use, e.g.,

    meson setup build --prefix=$PWD/build-install

The installation can be tested with

    meson test -C build

Note that the instructions above assume that all dependencies are found. For the C code there are only a few requirements besides meson itself, and that is basically BLAS and FFTW. For the rest see `l2cacc.yml` and `l2copenblas.yml` for two lists of dependencies, where the first makes use of the Accelerate framework and native compilers on a MacBook Pro M3. The `l2copenblas.yml` is more generic and pulls in everything from Conda, including compilers and OpenBlas. You can set up the environment using

    conda env --create -f l2copenblas.yml
    conda activate l2copenblas

You should then

    cd src
    ./build_meson.sh
    export PATH=$PWD/build-install/bin:$PATH
    export LD_LIBRARY_PATH=$PWD/build-install/lib:$LD_LIBRARY_PATH

and you should then be ready to run the l2c executable, for example

    l2c -N512 -d2

to run a L2C followed by a C2L, checking for accuracy.

# Codespace
Another simple way to test this code is to create a codespace. The l2copenblas.yml file in the root folder will then make sure that the codespace creates a conda environment with all necessary dependencies, including OpenBlas and FFTW, already installed. Just press the codespace button and wait awhile for the environment to build. Then enable the environment and run some tests or test the executable `l2c`

     bash # Create new terminal that is set up to run l2c
     cd src
     l2c -N1000 -d2 # runs a forward and backward transform and computes the error for an array of length 1000
     meson test -C build # runs all the tests
     cd results
     python table1.py # Create Table 1 in the paper
     python figure3.py # Create Figure 3 in the paper
     etc..
