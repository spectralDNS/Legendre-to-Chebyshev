l2c = executable('l2c',
                 'test.c',
                 include_directories : '..',
                 dependencies : [cblas, fftw, omp],
                 link_with : l2clib)
test('l2c', l2c)
