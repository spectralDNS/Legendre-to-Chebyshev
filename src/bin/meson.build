prefix = get_option('prefix')
lib_dir = prefix / get_option('libdir')
l2c = executable('l2c',
                 'l2c.c',
                 include_directories : '../C/',
                 dependencies : deps,
                 install : true,
                 #c_args : ['-march=native', '-Ofast'],
                 install_rpath : lib_dir,
                 link_with : l2clib)

test('test l2c forward-backward N1000', l2c, args : ['-N1000', '-d2'])
test('test l2c forward-backward N100', l2c, args : ['-N100', '-d2'])
test('test l2c direct', l2c, args : ['-N100', '-d3'])
test('test l2c forward 2d', l2c, args : ['-N400', '-d4'])
test('test l2c forward-backward 2d', l2c, args : ['-N400', '-d5'])
test('test l2c dct', l2c, args : ['-N400', '-d4'])