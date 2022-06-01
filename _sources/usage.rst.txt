Usage
=====
Instructions for installation and use. Issues may be reported at https://github.com/uiuc-ae-l3harris/APC-L3Harris.

Installation from source
------------------------
Prerequisites
_____________

* A C++ compiler that supports the C++11 standard as well as OpenMP
* Python 3.6+
* OpenMP runtime library

Preparing the source code
_________________________
Clone the APC github repository to your computer. The root of this directory will be referred to as 'APC'. Create a directory in 'APC' called 'extern'. The root folder structure should be similar to:
::

    APC
    ├── CMakeLists.txt
    ├── README.txt
    ├── UserGuide.docx
    ├── UserGuide.pdf
    ├── bin
    ├── docs_doxygen
    ├── docs_sphinx
    ├── extern
    ├── getspice.py
    ├── include
    ├── pyproject.toml
    ├── setup.py
    ├── src
    ├── srcMB
    ├── srcPy
    └── srcTest
  
    
Go to https://naif.jpl.nasa.gov/naif/toolkit_C.html and download the correct cspice library for your system. The file will be called either "cspice.tar.z" or "cspice.zip". Extract this file into the 'APC/extern' folder. The Folder structure should now look like:

::

    APC
    ├── CMakeLists.txt
    ├── README.txt
    ├── UserGuide.docx
    ├── UserGuide.pdf
    ├── bin
    ├── docs_doxygen
    ├── docs_sphinx
    ├── extern
    │   └── cspice
    │       ├── N0067
    │       ├── data
    │       ├── doc
    │       ├── etc
    │       ├── exe
    │       ├── include
    │       ├── lib
    │       ├── makeall.csh
    │       └── src
    ├── getspice.py
    ├── include
    ├── pyproject.toml
    ├── setup.py
    ├── src
    ├── srcMB
    ├── srcPy
    └── srcTest

Preparing the Python virtual environment
________________________________________

To prevent conflicting prerequisites, work in Python should be performed in a virtual environment. I recommend Anaconda for managing virtual environments. Alternatively, in Python 3.6+ the venv command can be used https://docs.python.org/3/library/venv.html.

Once your virtual environment has been setup and activated install PyBind11 with pip.

.. code-block:: console

   (.venv) APC $ pip install pybind11

.. code-block:: console

   (.venv) APC $ conda install -c conda-forge pybind11

If using pip or conda to install python packages is not an option then go to the 'Installing PyBind11 Manually'_.

.. _Building APC:

Building APC
____________

Move your consoles current working directory to "APC" and build the APC library.

.. code-block:: console
    
    (.venv) APC $ python setup.py build_ext --inplace

Once the build finishes, a python library similar to "APC.cpython-38-x86_64-linux-gnu.so" should be present in "APC". This library can be moved to your desired working directory and imported into python by

.. code-block:: python

    import APC

Building APC Alternative
________________________
Instead of building the library inplace, the APC library can be installed in the current virtual environment with pip. In "APC" run the following

.. code-block:: console

    (.venv) APC $ pip install .

APC will be installed as a package in the virtual environment and may be imported without the need to have the library in the working directory so long as the virtual environment is active.

To uninstall APC, run the following while the virtual environment is active

.. code-block:: console
    
    (.venv) $ pip uninstall APC

.. _Installing Pybind11 Manually:

Installing Pybind11 Manually
----------------------------
This method will work for compiling APC with a PyBind11 wrapper as long as you have access to Github. This section assumes you have followed the Installation from Source guide up to creating the python virtual environment.

In the root folder where the APC repo was cloned run this git command to switch to the includePybind11 branch.

.. code-block:: console
    
    (.venv) $ git checkout includePybind11

A few files will change and in the extern folder there will be a submodule folder called "PyBind11" that will be empty. To sync this folder with the main PyBind11 repo run

.. code-block:: console
    
    (.venv) $ git submodule update --init

The folder should populate with the PyBind11 repo. Go to 'Building APC'_ to finish compiling the APC library from source.

Usage examples
==============

In the "bin" folder within "APC" there are two python notebooks that make use of the APC library. Test.ipynb demonstrates APC running sequentially and in parallel with a few orbtial plots. Benchmark.ipynb demonstrates the performance gains from running multiple orbits in parallel.
Running these notebooks requires the following Python packages:

* Numpy
* Plotly

For APC to utilize CSPICE, two CSPICE kernels must be present in the working directory

1. An Earth, Moon, Sun ephemeris kernel file. de440.bsp from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/ is suitable.
2. A leapseconds kernel file. naif0012.tls from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/ is suitable.

The matrices folder must also be present in the working directory. These files contain data for the 2008 Earth gravity model and has the structure:

::

    matrices
    ├── A_matrices.bin
    ├── P1_matrices.bin
    ├── P2_matrices.bin
    ├── T1_matrices.bin
    ├── T2_matrices.bin
    └── Ta_matrices.bin