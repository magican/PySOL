# encoding: utf-8
"""
Setup script for cerbere package. See package documentation for information.

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
import os
import sys
from distutils.core import setup

from numpy import get_include

with open("CONTRIBUTORS.txt", "r") as f:
    authors_email = ', '.join(f.readlines())

setup(
    name = "cerbere",
    version = "0.1.1",
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
    install_requires=[ 'numpy>=1.7.1'
                     , 'scipy>=0.12.0'
                     , 'netCDF4>=1.0.4'
                     , 'pyhdf==0.8.3'
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
