prefix = get_option('prefix')
lib_dir = prefix / get_option('libdir')
l2cprec = executable('l2cprec',
                     'L2Cprec.cpp',
                     include_directories : '../../C',
                     link_with : l2clib,
                     install_rpath : lib_dir,
                     install : true,
                     dependencies : deps) # hack link to get include directories
test('test multiprecision boost',
     l2cprec, args : ['-N100', '-m0', '-a0'])
