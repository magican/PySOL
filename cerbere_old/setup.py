# encoding: utf-8
"""
Setup script for cerber package. See package documentation for information.

This script is modelled on the Django install script.

:copyright: Copyright 2013 Pelamis Scientific Software Ltd
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: David Poulter <david.poulter@pelamis.co.uk>
.. codeauthor:: David Poulter <david.poulter@pelamis.co.uk>
"""
import os
import sys
from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
from numpy import get_include

# The name an location of Cython modules
#EXT_MODULES = [Extension("cerbere.utils.fastregrid", ["cerbere/utils/fastregrid.pyx"]),
#	       Extension("cerbere.utils.inpaint", ["cerbere/utils/inpaint.pyx"])]

with open("CONTRIBUTORS.txt", "r") as f:
    authors_email = ', '.join(f.readlines())

setup(
    name = "cerbere",
    version = "0.1.0",
    author = "Felyx Project Collaborators: %s" % authors_email,
    author_email = "jean.francois.piolle@ifremer.fr",
    description = "",
    license = "LICENSE.txt",
    keywords = "",
    url = "http://hrdds.ifremer.fr",
    packages = [ 'cerbere'
               , 'cerbere.datamodel'
               , 'cerbere.geo'
               , 'cerbere.mapper'
               , 'cerbere.mapper.test'
               ],
    package_data = {'cerbere': [ 'geo/resources/*.*'
                               ]
                   },
    long_description="",
    install_requires=[ 'Cython>=0.19.1'
                     , 'numpy>=1.7.1'
                     , 'scipy>=0.12.0'
                     , 'netCDF4>=1.0.4'
                     , 'pyhdf>=0.8.3'
                     , 'pygrib>=1.9.6'
                     , 'Shapely==1.2.18'
                     , 'python-dateutil>=2.1'
                     , 'GDAL>=1.7.0'
                     , 'pyresample>=1.1.0'
                     ],
    classifiers=[
        "Development Status :: 1 - Pre-Alpha",
        "License :: LICENSE.txt",
        ],
    #cmdclass = {'build_ext': build_ext},
    #ext_modules = EXT_MODULES,
    include_dirs=[get_include()]
    #cmdclass={"install_data": post_install}
)
