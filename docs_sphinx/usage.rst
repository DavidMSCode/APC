Usage
=====
Instructions for installation and use are a work in progress. Issues may be reported at https://github.com/DavidMSCode/APC/issues.

Installation from source
------------------------
Clone the APC github repository to your computer. The root of this directory will be referred to as 'APC.' Create a directory in 'APC' called 'extern.' The root folder structure should be similar to:
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


Next install APC using pip in the root directory of the APC source code:

.. code-block:: console

   (.venv) APC $ pip install .

Then in a python file or notebook import the module:

.. code-block:: python
    
    import APC

