#!/bin/sh

#------------------------------------------------------------------------------
# Cerbere dependencies installation for Ubuntu 12.04
#------------------------------------------------------------------------------
here="$(pwd)"

sudo apt-get install gcc \
                     g++ \
                     gfortran \
                     libblas-dev \
                     libatlas-dev \
                     liblapack-dev \
                     libhdf5-serial-dev \
                     libgdal1-dev \
                     libgrib-api-dev \
                     python-dev \
                     python-pip

# Numpy installation fails when it is done during "install_requires" processing
# of setup.py but it works fine when done this way.
pip install numpy==1.8

# the pyhdf project changed the version number convention they use for the
# python packages (X.Y-Z to X.Y.Z), which prevents pip from detecting the latest
# version correctly => manual installation.
cd /tmp
wget "http://downloads.sourceforge.net/project/pysclint/pyhdf/0.8.3/pyhdf-0.8.3.tar.gz"
tar xvzf pyhdf-0.8.3.tar.gz
cd pyhdf-0.8.3
if [ -f /usr/lib/libmfhdfalt.so ]; then
    # Patch setup.py to use the -alt version of the libraries
    cp "${here}"/pyhdf_setup.py.patch .
    patch -p0 < pyhdf_setup.py.patch
fi
export INCLUDE_DIRS=/usr/include/hdf
export LIBRARY_DIRS=/usr/lib
export NOSZIP=1
python setup.py install

# Now we can install cerbere
