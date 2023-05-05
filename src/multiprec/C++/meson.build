l2cprec = executable('l2cprec',
                     'L2Cprec.cpp',
                     link_with : l2clib)
test('test multiprecision boost',
     l2cprec, args : ['-N100', '-m0'])
