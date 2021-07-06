import sys
import os
import Cython
import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize
from nvcc_integration import locate_cuda, cuda_get_config_vars, cuda_CCompiler
import distutils.ccompiler as dccm

Cython.Compiler.Options.embed = True

include_dirs = os.getenv('IDIR')
include_dirs = include_dirs.replace("-I", "").split()
include_dirs.append(numpy.get_include())

source_files = os.getenv('SOURCE_FILES')
source_files =source_files.replace('./src', '../src').split()
source_files.append("py_corona_sim.pyx")
# print(source_files)

libraries = os.getenv('LIBRARIES')
libraries = libraries.split()
libraries.append('-lgfortran')

compile_flags = os.getenv('COMPILE_FLAGS')
compile_flags = compile_flags.split()
print(compile_flags)

extension = Extension("py_corona_sim",
                      sources=source_files,
                      extra_objects=['../bin/ipbackgroundCFR_fun.o'],
                      language="c++",
                      extra_link_args=libraries,
                      extra_compile_args=compile_flags,
                      include_dirs=include_dirs)


# relies on command line -D RT_FLOAT to compile GPU code
if "-RT_FLOAT" in sys.argv:
    print("compiling with Real=float")
    cython_compile_time = {'RT_FLOAT': True}
    sys.argv.remove("-RT_FLOAT")

    for arg in os.getenv('CUDA_DLTO').split():
        extension.extra_link_args.append(arg)

    import distutils.sysconfig as dsc
    dsc.get_config_vars = cuda_get_config_vars
    dccm.CCompiler = cuda_CCompiler

else:
    print("compiling with Real=double")
    cython_compile_time = {'RT_FLOAT': False}

    
# monkey-patch for parallel compilation
def parallelCCompile(self, sources,
                     output_dir=None,
                     macros=None,
                     include_dirs=None,
                     debug=1,
                     extra_preargs=None,
                     extra_postargs=None,
                     depends=None):
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N=20 # number of parallel compilations
    import multiprocessing.pool
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile,objects))
    return objects


dccm.CCompiler.compile = parallelCCompile

setup(name="py_corona_sim",
      author="Mike Chaffin",
      version="0.0",
      ext_modules=cythonize([extension],
                            force=True,
                            compile_time_env=cython_compile_time),
      zip_safe=False
      )
