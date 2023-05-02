project('leg2cheb', 'c', 'cython',
        version : '1.0',
        default_options : [
          'buildtype=debugoptimized',
          'warning_level=3'])

cy = meson.get_compiler('cython')
cython = find_program(cy.cmd_array()[0])

cblas = dependency('openblas')
fftw = dependency('fftw3')
deps = [cblas, fftw]

subdir('C')
subdir('cython')
subdir('tests')
subdir('bin')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : l2clib,
                 version : '1.0',
                 name : 'leg2cheb',
                 description : 'A Library to transform between Legendre and Chebyshev coefficients')
