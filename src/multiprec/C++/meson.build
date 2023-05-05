add_global_arguments('-I/usr/share/miniconda/envs/l2c/include', language : 'cpp')
l2cprec = executable('l2cprec',
                     'L2Cprec.cpp')
test('test multiprecision boost',
     l2cprec, args : ['-N100', '-m0'])
