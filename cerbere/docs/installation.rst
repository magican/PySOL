====================
Cerbere installation
====================


Requirements
============

* BLAS

* ATLAS

* Lapack

* GRIB API

* HDF5

* HDF4

* Netcdf4

* GDAL

* Python

* Python pip

* gfortran

* gcc

* g++


Install these with your package manager (if a "dev" version of the package is available, use it).

Your Python package manager should handle most of the dependencies by himself, but there are some gotchas:

- Numpy is required by the setup.py scripts of other modules, so it must be installed beforehand (see https://github.com/pypa/pip/issues/25)

  .. sourcecode :: bash
  
    pip install numpy==1.8.0

- The versioning convention of the pyhdf package changed from major.minor-revision to major.minor.revision. This change prevents package managers to handle pyhdf correctly, you have to install it manually

  .. sourcecode :: bash
  
    cd /tmp
    wget "http://downloads.sourceforge.net/project/pysclint/pyhdf/0.8.3/pyhdf-0.8.3.tar.gz"
    tar xvzf pyhdf-0.8.3.tar.gz
    cd pyhdf-0.8.3
    
    # Change the following variables depending on your environment
    export INCLUDE_DIRS=/usr/include/hdf
    export LIBRARY_DIRS=/usr/lib
    export NOSZIP=1
    
    python setup.py install

  .. note::

    If the installation fails, check that the directory you specified in LIBRARY_DIRS contains libmfhdf.so and libdf.so.
    
    They may have been renamed, preventing the linker to find them.
    
    For example, on Ubuntu 12.04, you can find these libraries as "libmfhdfalt.so" and "libdfalt.so".
    
    To fix this, edit setup.py and set the correct names in the "libraries" variable (line 88).


On Ubuntu 12.04, you can use the install.sh script provided in the cerbere package to install dependencies.


Cerbere
=======

Just use your package manager to perform installation

  .. sourcecode :: bash
  
    pip install ./cerbere-0.0.1dev.tar.gz

  .. note::

    If pip complains about a missing C file, you should try to install / update
    distribute:

    .. sourcecode:: bash

      pip install -U distribute