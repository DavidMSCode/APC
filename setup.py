
from distutils.command.build_ext import build_ext
from glob import glob
from pathlib import Path
import platform
from setuptools import setup
import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "extern", "pybind11"))
from pybind11.setup_helpers import Pybind11Extension        #gets setup helpers from manually copied PyBind11 or from Pip installed
del sys.path[-1]

__version__ = "0.0.2"
cwd = os.path.curdir
cspice_lib_dir = "extern/cspice/lib/"
# Compiler args
cpp_extra_args = []
link_args = []
if platform.system() == 'Darwin':
    # Mac OS unique args
    cpp_extra_args.append('-mmacosx-version-min=10.15')
    cpp_extra_args.append("-Xpreprocessor") # Enable OpenMP support on Clang++ for newer versions of Mac OS
    cpp_extra_args.append("-fopenmp")   #enable openmp
    cpp_extra_args.append('-std=c++11') #c++11 standard
    link_args.append('-mmacosx-version-min=10.15')
    link_args.append('-lomp')  # use llvm OpenMP runtime
    cspice_lib_name = "cspice.a" # link against cspice library (external to normal lib locations)

elif platform.system() == "Windows":
    # Windows Args
    link_args.append("/NODEFAULTLIB:msvcrt.lib")  #ignore incompatible runtime libs
    link_args.append("/NODEFAULTLIB:libcmtd.lib")
    link_args.append("/NODEFAULTLIB:msvcrtd.lib")
    cspice_lib_name = "cspice.lib"
    cpp_extra_args.append("/std:c++14") #visual studio doesn't include c++11 standard so use c++14 instead  
    cpp_extra_args.append("/openmp") #enable openmp

elif platform.system() == "Linux":
    # Linux args
    cspice_lib_name = "cspice.a"# link against cspice library (external to normal lib locations)
    cpp_extra_args.append("-fopenmp") # Enable openmp
    link_args.append('-lgomp')  # use GNU OpenMP library
    cpp_extra_args.append('-std=c++11')  # use C++11 standard
else:
    raise Exception("Attempting to build APC python package with unsupported OS.")

"""Check if a static cspice library is located in extern or if a shared library should be built"""
if Path(os.path.join(cspice_lib_dir,cspice_lib_name)).is_file():
    print("User provided static CSPICE library")
else:
    if platform.system()=="Darwin" or platform.system()=="Linux":
        cspice_lib_name = "libcspice.so"
    else:
        cspice_lib_name = "libcspice.dll"
    #build library if get_cspice has not already built cspice library
    if not Path(os.path.join(cspice_lib_dir,cspice_lib_name)).is_file():
        print("Trying to build CSPICE library")
        from get_spice import main
        main()
    else:
        print("Shared CSPICE library already built")
        pass

cspice_lib_path = os.path.join(cspice_lib_dir,cspice_lib_name)


# Get all necessary C++ source files
SRCFOLDERS = ["src", "srcPy"]# src folder names
SRCFILES = []
for folder in SRCFOLDERS:
    SRCFILES = SRCFILES+glob(folder+"/*.cpp") # grab all .cpp files at each location
SRCFILES = sorted(SRCFILES) # Sort source files for consistency

"""Get header files"""
includes = []
includes.append("include") #APC header files location
includes.append("extern/cspice/include") #cspice library header files
if os.path.isdir("extern/pybind11/pybind11"):
    includes.append("extern/pybind11/include") #pybind11 header files (if included as external directory)
    print("Using PyBind11 library from extern/")

# Define the module
ext_modules = [
    Pybind11Extension(
        "APC",  # The module name
        SRCFILES,  # .cpp files to be compiled
        include_dirs=includes,  # header file locations for src files and external libraries
        extra_compile_args=cpp_extra_args,  # extra cpp compile args
        # link cspice library https://stackoverflow.com/questions/4597228/how-to-statically-link-a-library-when-compiling-a-python-module-extension
        extra_objects=[cspice_lib_path],
        extra_link_args=link_args  # Link to libs
    ),
]

# Setup module
setup(
    name="APC",
    version=__version__,
    cmdclass={"build_ext": build_ext},
    author="David Stanley",
    author_email="davidms4@illinois.edu",
    description="Test compilation",
    ext_modules=ext_modules,
    zip_safe=False,
    python_requires=">=3.6",
)
