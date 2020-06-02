from setuptools import setup, Extension
from Cython.Build import cythonize

# CPU version
extension = Extension("py_corona_sim",
                      ["py_corona_sim_cpu.pyx"],
                      extra_objects=['build/libobservation_fit.a'],
                      language="c++",
                      extra_compile_args=['-lm','-fopenmp','-O3','-march=native','-DNDEBUG','-fPIC'],
                      include_dirs=["../src/",
                                    "/home/mike/Documents/Utilities/boost_1_73_0/",
                                    "/home/mike/Documents/Utilities/eigen-3.3.7/"])

setup(name="py_corona_sim",
      author = "Mike Chaffin",
      version = "0.0",
      
      ext_modules = cythonize([extension]),

      zip_safe = False
      )
