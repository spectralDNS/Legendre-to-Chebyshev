add_global_arguments('-DTRI', language : 'c')
add_global_link_arguments('-lm', language : 'c')
install_headers('leg2cheb.h')
l2clib = library('leg2cheb',
                 'leg2cheb.c',
                 install : true,
                 version : '1.0',
                 dependencies : deps)
