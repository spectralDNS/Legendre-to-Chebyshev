l2cprec = executable('l2cprec',
                     'L2Cprec.cpp',
                     dependencies : cblas) # hack link to get include directories
test('test multiprecision boost',
     l2cprec, args : ['-N100', '-m0'])
