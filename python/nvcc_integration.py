# from future.utils import iteritems
import os
from os.path import join as pjoin
from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


def find_in_path(name, path):
    """Find a file in a search path"""

    # Adapted fom http://code.activestate.com/recipes/52224
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system
    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.
    Starts by looking for the CUDA_HOME env variable. If not found,
    everything is based on finding 'nvcc' in the PATH.
    """

    # First check if the CUDA_HOME env variable is in use
    if 'CUDA_HOME' in os.environ:
        home = os.environ['CUDA_HOME']
        nvcc = pjoin(home, 'bin', 'nvcc')
    else:
        # Otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            raise EnvironmentError('The nvcc binary could not be '
                'located in your $PATH. Either add it to your path, '
                'or set $CUDA_HOME')
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home': home, 'nvcc': nvcc,
                  'include': pjoin(home, 'include'),
                  'lib64': pjoin(home, 'lib64')}
    for k, v in iter(cudaconfig.items()):
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be '
                                   'located in %s' % (k, v))

    return cudaconfig

# The following functions pull compiler flags out of 
from distutils.sysconfig import get_config_vars as default_get_config_vars

wrap_nvcc_flags = ['-Werror=format-security',
                   '-Wno-unused-result',
                   '-Wsign-compare',
                   '-Wformat',
                   '-D_FORTIFY_SOURCE=2',
                   '-fPIC',
                   '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']

remove_nvcc_flags = ['-O2',
                     '-Wall',
                     '-fwrapv',
                     '-Wdate-time',
                     '-fstack-protector-strong',
                     '-pthread']

def remove_compiler_warning_flags(x):
    if isinstance(x, str):
        for f in remove_nvcc_flags:
            x = x.replace(f, '')
        for f in wrap_nvcc_flags:
            x = x.replace(f, '-Xcompiler '+f)
        x = x.replace("-Wl,", "-Xlinker ")
    return x

def cuda_get_config_vars(*args):
    # print(args)
    result = default_get_config_vars(*args)
    # sometimes result is a list and sometimes a dict:
    if isinstance(result, list):
        return [remove_compiler_warning_flags(x) for x in result]
    elif isinstance(result, dict):
        return {k: remove_compiler_warning_flags(x)
                for k, x in result.items()}
    else:
        raise Exception("cannot handle type"+type(result))
