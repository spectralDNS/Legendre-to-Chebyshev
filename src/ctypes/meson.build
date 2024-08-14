# needs LD_LIBRARY_PATH and/or DYLD_LIBRARY_PATH
py = import('python').find_installation()
py.install_sources(['l2ctypes.py'])

test('test ctypes',
     py,
     args: [
       'l2ctypes.py',
     ],
     depends : l2clib,
     workdir: meson.current_source_dir())