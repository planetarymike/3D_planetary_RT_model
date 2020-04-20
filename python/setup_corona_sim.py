# setup_multi.py --- setup routine for wrapping C++ H corona code

from setuptools import setup, Extension
from Cython.Build import cythonize

extension = Extension("py_corona_sim",
                      ["py_corona_sim.pyx"],      # our Cython source
                      language="c++",             # generate C++ code
                      extra_compile_args=["-std=c++17","-O3",
                                          "-ftree-vectorize","-march=native",
                                          "-ffast-math","-funsafe-math-optimizations","-mfpmath=sse",
                                          "-DNDEBUG"],
                      extra_link_args=["-lm","-fPIC","-fopenmp"],
                      include_dirs=["../bin/",
                                    "/home/mike/Documents/Utilities/boost_1_72_0/",
                                    "/home/mike/Documents/Utilities/eigen-3.3.7/"])

#    extra_objects=["/required_procedures/ipbackgroundCFR_fun.o"]

setup(name="py_corona_sim",
      ext_modules = cythonize(extension),
      )


