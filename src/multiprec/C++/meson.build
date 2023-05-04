boost = cpp.find_library('boost_system')
l2cprec = executable('l2cprec',
                     'L2Cprec.cpp',
                     include_directories : '.',
                     dependencies : boost)
test('test multiprecision boost', l2cprec, args : ['-N100', '-m0'])
