install_headers('leg2cheb.h')
l2clib = library('leg2cheb',
                 'leg2cheb.c',
                 install : true,
                 c_args : ['-march=native', '-shared'],
                 version : '1.0',
                 dependencies : deps)
