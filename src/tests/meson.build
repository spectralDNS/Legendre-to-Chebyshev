l2ctest = executable('l2ctest',
                     'test.c',
                     include_directories : '../C',
                     dependencies : deps,
                     link_with : l2clib)
test('test libleg2cheb', l2ctest)

py = import('python').find_installation()
test('test cython',
     py,
     args: [
       'test.py',
     ],
     workdir: meson.current_source_dir())
