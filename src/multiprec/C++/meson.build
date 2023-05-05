l2cprec = executable('l2cprec',
                     'L2Cprec.cpp',
                     dependencies : cblas)
test('test multiprecision boost',
     l2cprec, args : ['-N100', '-m0'])
