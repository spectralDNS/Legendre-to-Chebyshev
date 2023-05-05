l2cprec = executable('l2cprec',
                     'L2Cprec.cpp',
                     include_directories : '/usr/share/miniconda/envs/l2c/include')
test('test multiprecision boost',
     l2cprec, args : ['-N100', '-m0'])
