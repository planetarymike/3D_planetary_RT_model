# setup_multi.py --- setup routine for wrapping C++ H corona code

from setuptools import setup, Extension
from Cython.Build import cythonize
#from nvcc_integration import *

# extension = Extension("py_corona_sim",
#                       ["py_corona_sim.pyx"],
#                       library_dirs = [CUDA['lib64']],
#                       libraries = ['cudart'],
#                       language="c++",
#                       runtime_library_dirs = [CUDA['lib64']],
#                       extra_compile_args=['-std=c++17','-lm','-fopenmp','-O3','-march=native','-DNDEBUG','-fPIC'],
#                       extra_objects=["observation_fit.o","observation_fit_gpu_device.o"],
#                       include_dirs=[CUDA['include']]
#                       )

extension = Extension("py_corona_sim",
                      ["py_corona_sim.pyx"], 
                      language="c++",
                      extra_compile_args=["-std=c++17","-O3","-march=native","-DNDEBUG","-fPIC"],
                      extra_link_args=["-lm","-fopenmp"],
                      include_dirs=["../bin/",
                                    "/home/mike/Documents/Utilities/boost_1_72_0/",
                                    "/home/mike/Documents/Utilities/eigen-3.3.7/"])


setup(name="py_corona_sim",
      author = "Mike Chaffin",
      version = "0.0",
      
      ext_modules = cythonize([extension]),

      zip_safe = False
      )
