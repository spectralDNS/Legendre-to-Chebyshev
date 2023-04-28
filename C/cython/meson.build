py = import('python').find_installation()
dep_py = py.dependency()
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

py.extension_module(
    'l2c',
    'l2c.pyx',
    install : true,
    include_directories: incdir_numpy,
    link_args : ['-lm'],
    dependencies : [dep_py,fftw,cblas],
    link_with : l2clib
)
