from distutils.command.build_ext import build_ext
from glob import glob
import platform
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
import os


__version__ = "0.0.1"
cpp_extra_args = []
link_args = []
cspice_lib = []

#Compiler args
if platform.system() == 'Darwin':
   cpp_extra_args.append('-mmacosx-version-min=10.9')

if os.name == "nt":
    link_args.append("/NODEFAULTLIB:library")
    cspice_lib = "extern/cspice/lib/cspice.lib"
else:
    cspice_lib = "extern/cspice/lib/cspice.a"
    cpp_extra_args.append('-std=c++11')  

#Get all necessary C++ source files
SRCFOLDERS = ["src","srcPy"]
SRCFILES=[]
for folder in SRCFOLDERS:
    SRCFILES=SRCFILES+glob(folder+"/*.cpp")
SRCFILES = sorted(SRCFILES)                                         # Sort source files for consistency

#Define external module
ext_modules = [
    Pybind11Extension(
        "APC",
        SRCFILES,                                                   #.cpp files to be compiled
        include_dirs=["include","extern/cspice/include"],           #header file locations for src files and external libraries
        extra_compile_args=cpp_extra_args,                          #extra cpp compile args
        extra_objects=[cspice_lib],                                 #link cspice library https://stackoverflow.com/questions/4597228/how-to-statically-link-a-library-when-compiling-a-python-module-extension
        extra_link_args=link_args
    ),
]

#Setup module
setup(
    name = "APC",
    version = __version__,
    cmdclass={"build_ext": build_ext},
    author="David Stanley",
    author_email = "davidms4@illinois.edu",
    description = "Test compilation",
    ext_modules = ext_modules,
    zip_safe = False,
    python_requires = ">=3.6",
)