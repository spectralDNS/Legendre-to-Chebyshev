py = import('python').find_installation()
#py.install_sources(['Leg2Cheb.py'])

test('test pure python',
     py,
     args: [
       'Leg2Cheb.py',
     ],
     workdir: meson.current_source_dir())

test('test cython',
     py,
     args: [
       'test.py',
     ],
     workdir: meson.current_source_dir())
