project('leg2cheb', 'c', 'cython', 'cpp',
        version : '1.0',
        default_options : [
          'buildtype=debugoptimized',
          'warning_level=3'])

cy = meson.get_compiler('cython')
cpp = meson.get_compiler('cpp')
cython = find_program(cy.cmd_array()[0])

cblas = dependency('openblas')
fftw = dependency('fftw3')
deps = [cblas, fftw]
add_global_arguments('-DTRI', language : 'c')
add_global_link_arguments('-lm', language : 'c')
add_global_arguments('-std=c++11', language : 'cpp')

subdir('C')
subdir('cython')
subdir('tests')
subdir('bin')
subdir('multiprec/python')
subdir('multiprec/C++')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : l2clib,
                 version : '1.0',
                 name : 'leg2cheb',
                 description : 'A Library to transform between Legendre and Chebyshev coefficients')
