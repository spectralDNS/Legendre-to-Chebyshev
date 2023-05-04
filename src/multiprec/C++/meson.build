boost = dependency('boost')
l2cprec = executable('l2cprec',
                     'L2Cprec.cpp',
                     dependencies : [boost])
test('test multiprecision boost',
     l2cprec, args : ['-N100', '-m0'])
