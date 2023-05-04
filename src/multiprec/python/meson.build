py = import('python').find_installation()
#py.install_sources(['L2Cprecision.py'])

test('test multiprecision python',
     py,
     args: [
       'L2Cprecision.py',
     ],
     workdir: meson.current_source_dir())
