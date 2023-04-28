project('leg2cheb', 'c', 'cython',
        version : '1.0',
        default_options : [
          'buildtype=debugoptimized',
          'warning_level=3'])

cy = meson.get_compiler('cython')
cython = find_program(cy.cmd_array()[0])

cblas = dependency('openblas')
fftw = dependency('fftw3')
omp = dependency('openmp')
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)
deps = [cblas, fftw, omp, m_dep]
add_global_arguments('-DTRI', language : 'c')
add_global_arguments('-fno-builtin', language : 'cython')

install_headers('leg2cheb.h')

l2clib = library('leg2cheb',
                 'leg2cheb.c',
                 install : true,
                 dependencies : deps)

subdir('tests')
subdir('cython')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : l2clib,
                 version : '1.0',
                 name : 'leg2cheb',
                 description : 'A Library to transform between Legendre and Chebyshev coefficients')
