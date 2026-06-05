import sys
import os
import Cython
import numpy

from setuptools import setup, Extension
from Cython.Build import cythonize

if sys.version_info[0] == 3 and sys.version_info[1] < 12:
    import distutils.ccompiler as dccm
    import distutils.sysconfig as dsc
    from distutils.ccompiler import CCompiler as default_CCompiler
    from distutils.sysconfig import get_config_vars as default_get_config_vars
    from distutils.sysconfig import customize_compiler as default_customize_compiler
elif sys.version_info[0] == 3 and sys.version_info[1] >= 12:
    # TODO: import from Cython or setuptools instead
    import distutils.sysconfig as dsc
    from distutils.sysconfig import get_config_vars as default_get_config_vars
    from distutils.sysconfig import customize_compiler as default_customize_compiler

Cython.Compiler.Options.embed = True

# read environment variables set by makefile to determine configuration
include_dirs = os.getenv('IDIR')
if "-isystem" in include_dirs:
    sysinclude_dirs = [d.split()[0] for d in include_dirs.split("-isystem")[1:] if d!='']
else:
    sysinclude_dirs = []
if "-I" in include_dirs:
    include_dirs = [d.split()[0] for d in include_dirs.split("-I")[1:] if d!='']
else:
    include_dirs = []
sysinclude_dirs.append(numpy.get_include())
print(f"{sysinclude_dirs = }")
print(f"{include_dirs = }")

source_files = os.getenv('SOURCE_FILES')
source_files = source_files.replace('./src', '../src').split()
source_files.append("py_corona_sim.pyx")
# print(source_files)

libraries = os.getenv('LIBRARIES')
libraries = libraries.split()
libraries.append('-lgfortran')
print(f"{libraries = }")

compile_flags = os.getenv('COMPILE_FLAGS')
compile_flags = compile_flags.split()
for d in sysinclude_dirs:
    if os.environ['CC'] == 'nvcc':
        compile_flags.append("-isystem="+d)
    else:
        compile_flags.append("-isystem"+d)
print(f"{compile_flags = }")

# set name of extension based on selected compiler
extension_name = 'py_corona_sim'
if os.environ['CC'] == 'nvcc':
    extension_name += '_gpu'
    if os.getenv('NVCC_CCBIN') != "":
        compile_flags.append(f"-ccbin={os.getenv('NVCC_CCBIN')}")
else:
    extension_name += '_cpu'
# print(f'\n\n{extension_name = }\n\n')

# set up extension defaults
extension = Extension(extension_name,
                      sources=source_files,
                      extra_objects=['../bin/ipbackgroundCFR_fun.o'],
                      language="c++",
                      extra_link_args=libraries,
                      extra_compile_args=compile_flags,
                      include_dirs=include_dirs)


# The following functions help to overwrite setuptools._distutils defaults
nvcc_flags_to_wrap = ['-fno-strict-overflow',
                      '-fno-omit-frame-pointer',
                      '-mno-omit-leaf-frame-pointer',
                      '-fstack-clash-protection',
                      '-fcf-protection',
                      '-fPIC',
                      '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION']

compiler_flags_to_remove = ['-O2',
                            '-Wall',
                            '-Werror=format-security',
                            '-Wno-unused-result',
                            '-Wsign-compare',
                            '-Wformat',
                            '-D_FORTIFY_SOURCE=2',
                            '-fwrapv',
                            '-Wdate-time',
                            '-fstack-protector-strong',
                            '-pthread',
                            '-Wl,-O1',
                            '-Wl,-Bsymbolic-functions',
                            '-Wl,-z,relro',
                            '-fwrapv',
                            '-fstack-protector-strong',
                            '-Werror=format-security',
                            '-Wstrict-prototypes',
                            '-Wunreachable-code']


def remove_compiler_flags(x):
    if isinstance(x, str):
        for f in compiler_flags_to_remove:
            x = x.replace(f, '')
        if os.environ['CC'] == 'nvcc':
            for f in nvcc_flags_to_wrap:
                x = x.replace(f, '-Xcompiler '+f)
        x = " ".join([el for el in x.split() if (("-Wl" not in el)
                                                 and (el != '-g'))])
        if os.environ['CC'] == 'nvcc':
            x.replace("-c", "-dc")
    return x


def my_get_config_vars(*args):
    # print(args)
    result = default_get_config_vars(*args)

    # sometimes result is a list and sometimes a dict:
    if isinstance(result, list):
        return [remove_compiler_flags(x) for x in result]

    if isinstance(result, dict):
        return {k: remove_compiler_flags(x)
                for k, x in result.items()}

    raise Exception("cannot handle type "+type(result))


# override default distutils compiler flag generator
dsc.get_config_vars = my_get_config_vars


def nvcc_get_cc_args(self, pp_opts, debug, before):
    # works for unixccompiler, cygwinccompiler
    cc_args = pp_opts + ['-dc']
    if debug:
        cc_args[:0] = ['-g']
    if before:
        cc_args[:0] = before
    return cc_args

# set GPU/CPU options based on presence of command line -D RT_FLOAT
if "-RT_FLOAT" in sys.argv:
    print("compiling with Real=float")
    cython_compile_time = {'RT_FLOAT': True}
    sys.argv.remove("-RT_FLOAT")

    for arg in os.getenv('CUDA_DLTO').split():
        extension.extra_link_args.append(arg)
else:
    print("compiling with Real=double")
    cython_compile_time = {'RT_FLOAT': False}


cython_compile_time['CPP_GIT_HASH'] = os.environ['CPP_GIT_HASH']


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
    macros, objects, extra_postargs, pp_opts, build = \
        self._setup_compile(output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    # parallel code
    N = 20  # number of parallel compilations
    import multiprocessing.pool
    def _single_compile(obj):
        try: src, ext = build[obj]
        except KeyError: return
        replace_extra_postargs = extra_postargs.copy()
        if 'py_corona_sim.cpp' in src and os.environ['CC']=='nvcc':
            # for the cython-generated cpp file only, suppress nvcc
            # warnings 177 and 550 that result from automatic cython
            # code generation
            replace_extra_postargs += ['--diag-suppress', '177,550']

        self._compile(obj, src, ext, cc_args, replace_extra_postargs, pp_opts)

    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(N).imap(_single_compile, objects))
    return objects


# overwrite link step for CUDA link time optimization
def override_link_shared_object(self,
                                objects,
                                output_filename,
                                output_dir=None,
                                libraries=None,
                                library_dirs=None,
                                runtime_library_dirs=None,
                                export_symbols=None,
                                debug=0,
                                extra_preargs=None,
                                extra_postargs=None,
                                build_temp=None,
                                target_lang=None):

    replace_extra_postargs = extra_postargs.copy()
    if os.environ['CC'] == 'nvcc':
        # CUDA wants different architectecture specifications when
        # compiling the individual object files and when linking all of
        # the object files together into a single binary, so we need to
        # overwrite some of the supplied arguments here

        # print('\nBEFORE SUBSTITUTION:')
        # print(f'{replace_extra_postargs = }')
        replace_extra_postargs = [arg.replace('code=lto_', 'code=sm_')
                                  for arg in replace_extra_postargs]
        replace_extra_postargs += ['-dlto']
        # print('\nAFTER SUBSTITUTION:')
        # print(f'{replace_extra_postargs = }\n')

    self.link(self.SHARED_OBJECT,
              objects,
              output_filename,
              output_dir,
              libraries,
              library_dirs,
              runtime_library_dirs,
              export_symbols,
              debug,
              extra_preargs,
              replace_extra_postargs,
              build_temp,
              target_lang)


def my_customize_compiler(compiler):
    default_customize_compiler(compiler)
    compiler.compile = parallelCCompile.__get__(compiler)
    compiler.link_shared_object = override_link_shared_object.__get__(compiler)
    if os.environ['CC'] == 'nvcc':
        compiler._get_cc_args = nvcc_get_cc_args.__get__(compiler)


if sys.version_info[0] == 3 and sys.version_info[1] < 12:
    # needed for EMUS python 3.9 environment
    dccm.CCompiler.compile = parallelCCompile
    dccm.CCompiler.link_shared_object = override_link_shared_object
    if os.environ['CC'] == 'nvcc':
        dccm.CCompiler._get_cc_args = nvcc_get_cc_args
elif sys.version_info[0] == 3 and sys.version_info[1] >= 12:
    dsc.customize_compiler = my_customize_compiler

setup(name=extension_name,
      author="Mike Chaffin",
      version="0.0",
      ext_modules=cythonize([extension],
                            force=True,
                            compile_time_env=cython_compile_time,
                            gdb_debug=True),
      zip_safe=False)
