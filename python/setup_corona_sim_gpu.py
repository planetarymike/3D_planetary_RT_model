from setuptools import setup, Extension
from Cython.Build import cythonize
from nvcc_integration import *

# GPU version
extension = Extension("py_corona_sim",
                      ["py_corona_sim_gpu.pyx"],
                      extra_objects=['build/libobservation_fit.a'],
                      library_dirs = [CUDA['lib64']],
                      libraries = ['cudart'],
                      language="c++",
                      runtime_library_dirs = [CUDA['lib64']],
                      extra_compile_args=['-lm','-fopenmp','-O3','-march=native','-DNDEBUG','-fPIC','-D RT_FLOAT'],
                      include_dirs=["../src/",
                                    "/home/mike/Documents/Utilities/boost_1_73_0/",
                                    "/home/mike/Documents/Utilities/eigen-3.3.7/",
                                    CUDA['include']])

setup(name="py_corona_sim",
      author = "Mike Chaffin",
      version = "0.0",
      
      ext_modules = cythonize([extension]),

      zip_safe = False
      )
