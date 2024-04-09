project('leg2cheb', 'c', 'cython', 'cpp',
        version : '1.0',
        default_options : [
          'buildtype=release',
          'warning_level=1',
          'optimization=3'])

cy = meson.get_compiler('cython')
cpp = meson.get_compiler('cpp')
cython = find_program(cy.cmd_array()[0])

cblas = dependency('openblas')
fftw = dependency('fftw3')
omp = dependency('openmp', required : false)
deps = [cblas, fftw, omp]
#add_global_arguments(['-march=native'], language : 'c')
add_global_arguments(['-fomit-frame-pointer'], language : 'c')
add_global_link_arguments(['-lm'], language : 'c')
#add_global_link_arguments(['-framework', 'Accelerate'], language : 'c')
add_global_arguments('-std=c++14', language : 'cpp')
if omp.found()
  add_global_arguments(['-DOMP'], language : 'c')
endif

#subdir('simd')
subdir('C')
subdir('cython')
subdir('ctypes')
subdir('bin')
subdir('python')
subdir('multiprec/python')
subdir('multiprec/C++')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : l2clib,
                 version : '1.0',
                 name : 'leg2cheb',
                 description : 'A Library to transform between Legendre and Chebyshev coefficients')
