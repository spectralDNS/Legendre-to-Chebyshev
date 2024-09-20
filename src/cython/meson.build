py = import('python').find_installation()
dep_py = py.dependency()
incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

pyversion = run_command(py,
  ['-c', 'import sys; print("python"+f"{sys.version_info.major}.{sys.version_info.minor}")'],
  check : true
).stdout().strip()

py.extension_module(
    'l2c',
    'l2c.pyx',
    install : true,
    install_dir : get_option('prefix') / 'lib' / pyversion / 'site-packages',
    c_args : ['-march=native'],
    include_directories: incdir_numpy,
    dependencies : [dep_py, deps],
    link_with : l2clib
)
