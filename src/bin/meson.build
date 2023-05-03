prefix = get_option('prefix')
lib_dir = prefix / get_option('libdir')
l2c = executable('l2c',
                 'l2c.c',
                 include_directories : '../C/',
                 dependencies : deps,
                 install : true,
                 install_rpath : lib_dir,
                 link_with : l2clib)

test('test l2c', l2c, args : ['-N1000', '-d2'])
