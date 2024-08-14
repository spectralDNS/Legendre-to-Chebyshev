prefix = get_option('prefix')
lib_dir = prefix / get_option('libdir')
install_headers('leg2cheb.h', install_dir : prefix / 'include')
l2clib = library('leg2cheb',
                 ['leg2cheb.c', 'dct.c', 'lambda.c'],
                 install : true,
                 install_dir : prefix / 'lib',
                 #c_args : ['-march=native', '-Ofast'],
                 version : '1.0',
                 dependencies : deps)
