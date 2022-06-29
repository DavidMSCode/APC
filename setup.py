from distutils.command.build_ext import build_ext
from glob import glob
import platform
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension
import os


__version__ = "0.0.2"

#Compiler args
cpp_extra_args = []
link_args = []
cspice_lib = []
if platform.system() == 'Darwin':
    #Mac OS unique args
   cpp_extra_args.append('-mmacosx-version-min=10.15')
   cpp_extra_args.append("-Xpreprocessor")                          #Enable OpenMP support on Clang++ for newer versions of Mac OS
   cpp_extra_args.append("-fopenmp")
   cpp_extra_args.append('-std=c++11')
   link_args.append('-mmacosx-version-min=10.15') 
   link_args.append('-lomp')                                        #llvm OpenMP   
   cspice_lib = "extern/cspice/lib/cspice.a"                       #link against cspice library (external to normal lib locations)    
                
elif platform.system() == "Windows":
    #Windows Args
    link_args.append("/NODEFAULTLIB:library")
    cspice_lib = "extern/cspice/lib/cspice.lib" 
    cpp_extra_args.append("/std:c++11")                    

elif platform.system()=="Linux":
    #Linux args                                                               
    cspice_lib = "extern/cspice/lib/cspice.a"                       #link against cspice library (external to normal lib locations)
    cpp_extra_args.append("-fopenmp")                               #Enable multithreading through pragma statements
    link_args.append('-lgomp')                                       #GNU OpenMP library 
    cpp_extra_args.append('-std=c++14')                             #use C++11 standard
    

#Get all necessary C++ source files
SRCFOLDERS = ["src","srcPy"]                                        # src folder names
SRCFILES=[]
for folder in SRCFOLDERS:
    SRCFILES=SRCFILES+glob(folder+"/*.cpp")                         # grab all .cpp files at each location      
SRCFILES = sorted(SRCFILES)                                         # Sort source files for consistency

#Define the module
ext_modules = [
    Pybind11Extension(
        "APC",                                                      #The module name
        SRCFILES,                                                   #.cpp files to be compiled
        include_dirs=["include","extern/cspice/include"],           #header file locations for src files and external libraries
        extra_compile_args=cpp_extra_args,                          #extra cpp compile args
        extra_objects=[cspice_lib],                                 #link cspice library https://stackoverflow.com/questions/4597228/how-to-statically-link-a-library-when-compiling-a-python-module-extension
        extra_link_args=link_args                                   #Link to libs
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