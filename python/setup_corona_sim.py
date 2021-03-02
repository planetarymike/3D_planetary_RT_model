import sys
from setuptools import setup, Extension
from Cython.Build import cythonize

import numpy

from nvcc_integration import locate_cuda

boost_dir = "/home/mike/Documents/Utilities/boost_1_73_0/"
eigen_dir = "/home/mike/Documents/Utilities/eigen-3.3.7/"

extension = Extension("py_corona_sim",
                      ["py_corona_sim.pyx"],
                      extra_objects=['build/libobservation_fit.a',
                                     '../bin/ipbackgroundCFR_fun.o'],
                      language="c++",
                      extra_link_args=['-lgfortran', '-lgfortran', '-lm'],
                      extra_compile_args=['-fopenmp', '-O3', '-march=native',
                                          '-DNDEBUG', '-fPIC'],
                      include_dirs=[numpy.get_include(),
                                    "../src/",
                                    boost_dir,
                                    eigen_dir])

# relies on command line -D RT_FLOAT to compile GPU code
if "-RT_FLOAT" in sys.argv:
    print("compiling with Real=float")

    CUDA = locate_cuda()

    extension.extra_compile_args.append('-D RT_FLOAT')
    extension.libraries.append('cudart')
    extension.libraries.append('cudadevrt')
    extension.libraries.append('cusolver')
    extension.libraries.append('cublas')
    extension.library_dirs.append(CUDA['lib64'])
    extension.runtime_library_dirs.append(CUDA['lib64'])
    extension.include_dirs.append(CUDA['include'])

    cython_compile_time = {'RT_FLOAT': True}

    sys.argv.remove("-RT_FLOAT")
else:
    print("compiling with Real=double")
    cython_compile_time = {'RT_FLOAT': False}

setup(name="py_corona_sim",
      author="Mike Chaffin",
      version="0.0",
      ext_modules=cythonize([extension],
                            compile_time_env=cython_compile_time),
      zip_safe=False
      )
