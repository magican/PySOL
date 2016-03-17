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
try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess

from numpy import get_include

gdal_deps = 'gdal>=1.7.0'
try:
    output = subprocess.check_output(['gdal-config', '--version'])
    major, minor, _ = output.split('.')
    gdal_deps = 'gdal>={}.{},<{}.{}'.format(major, minor, major, (int(minor) + 1))
except subprocess.CalledProcessError:
    # Silenced
    pass

with open("CONTRIBUTORS.txt", "r") as f:
    authors_email = ', '.join(f.readlines())

# Update version when there is a new git commit
major_minor_version = '0.1'
package_dir = os.path.dirname(__file__)
version_path = os.path.join(package_dir, 'VERSION.txt')
if os.path.exists('.git') and os.path.isdir('.git'):
    commits = subprocess.check_output([ '/usr/bin/git'
                                      , 'rev-list'
                                      , 'HEAD'
                                      , '--count']).decode('utf-8').strip()
    with open(version_path, 'w') as f:
        f.write('{}.{}\n'.format(major_minor_version, commits))

with open(version_path, 'r') as f:
    version = f.read()

setup(
    name = "cerbere",
    version = version,
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
                     , 'python-hdf4>=0.9'
                     , 'pygrib>=1.9.6'
                     , 'Shapely==1.2.18'
                     , 'python-dateutil>=2.1'
                     , gdal_deps
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
