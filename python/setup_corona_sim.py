import sys
import os
import Cython
import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize
from nvcc_integration import locate_cuda, cuda_get_config_vars

include_dirs = os.getenv('IDIR')
include_dirs = include_dirs.replace("-I", "").split()
include_dirs.append(numpy.get_include())

extra_objects = os.getenv('COMPILED_OBJECTS')
extra_objects = extra_objects.replace('./bin', '../bin').split()
extra_objects.append('../bin/ipbackgroundCFR_fun.o')

extension = Extension("py_corona_sim",
                      ["py_corona_sim.pyx"],
                      extra_objects=extra_objects,
                      # extra_objects=['build/libobservation_fit.a',
                      #                '../bin/ipbackgroundCFR_fun.o'],
                      language="c++",
                      extra_link_args=['-lgfortran', '-lgfortran', '-lm'],
                      extra_compile_args=['-O3', '-DNDEBUG'],
                      include_dirs=include_dirs)


# relies on command line -D RT_FLOAT to compile GPU code
if "-RT_FLOAT" in sys.argv:
    print("compiling with Real=float")

    CUDA = locate_cuda()

    extension.extra_compile_args.append('-D RT_FLOAT')
    for arg in os.getenv('CUDA_DLTO').split():
        extension.extra_compile_args.append(arg)
        extension.extra_link_args.append(arg)
    # extension.extra_link_args.append('--verbose')
    extension.libraries.append('cudart')
    extension.libraries.append('cudadevrt')
    extension.libraries.append('cusolver')
    extension.libraries.append('cublas')
    extension.library_dirs.append(CUDA['lib64'])
    #extension.runtime_library_dirs.append(CUDA['lib64'])
    extension.include_dirs.append(CUDA['include'])

    cython_compile_time = {'RT_FLOAT': True}

    sys.argv.remove("-RT_FLOAT")

    import distutils.sysconfig as dsc
    dsc.get_config_vars = cuda_get_config_vars

else:
    print("compiling with Real=double")

    extension.extra_compile_args.append('-fopenmp')
    extension.extra_compile_args.append('-march=native')

    cython_compile_time = {'RT_FLOAT': False}

setup(name="py_corona_sim",
      author="Mike Chaffin",
      version="0.0",
      ext_modules=cythonize([extension],
                            force=True,
                            compile_time_env=cython_compile_time,
                            compiler_directives={'embedsignature': True}),
      zip_safe=False
      )
