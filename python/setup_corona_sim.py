# setup_multi.py --- setup routine for wrapping C++ H corona code

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(Extension(
    "py_corona_sim",
    ["py_corona_sim.pyx"],      # our Cython source
    language="c++",             # generate C++ code
    extra_compile_args=["-std=c++17","-O2","-ftree-vectorize","-march=native"],
    extra_link_args=["-lm","-fPIC"],
#    define_macros=[('SRCFNSLOC','"./source_functions/"')],
    include_dirs=["../bin/",
                  "/home/mike/Documents/Utilities/boost_1_72_0/",
                  "/home/mike/Documents/Utilities/eigen-3.3.7/"]))
#    extra_objects=["/required_procedures/ipbackgroundCFR_fun.o"]))
)
