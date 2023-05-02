l2ctest = executable('l2ctest',
                     'test.c',
                     include_directories : '../C',
                     dependencies : deps,
                     link_with : l2clib)
test('test lib', l2ctest)

py = import('python').find_installation()
test('pytests',
     py,
     args: [
       '-m', 'pytest',
       'test.py',
     ],
     workdir: meson.current_source_dir())
