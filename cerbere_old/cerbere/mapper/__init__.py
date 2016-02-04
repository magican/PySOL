# encoding: utf-8
__all__ = ['abstractmapper', 'ghrsstncfile', 'gribfile', 'hdffile',
    'ncfile', 'qscathdffile', 'saralncfile', 'sopranoncfile', 'urlseries',
    'ascatifrncfile']

from pkg_resources import WorkingSet, DistributionNotFound
import logging

# Import what mappers are avaialble.
packages = WorkingSet() 

try:
    packages.require('netCDF4')
    # from . import ghrsstncfile
    # from . import ncfile
    # from . import saralncfile
    # from . import sopranoncfile
    # from . import ascatifrncfile
except DistributionNotFound, exception:
    logging.warning('Python netCDF4 package is required for netCDF. ' 
        'No netCDF compatable modules are avaialble on this cerbere instance')
except ImportError:
    logging.exception('A netCDF mapper failed to load, and is unavailable:')
    
try:
    packages.require('pyhdf')
    # from . import hdffile
    # from . import qscathdffile
except DistributionNotFound, exception:
    logging.warning('Python pyhdf package is required for HDF. ' 
        'No HDF compatable modules are avaialble on this cerbere instance')
except ImportError:
    logging.exception('An HDF mapper failed to load, and is unavailable:')    
    
try:
    packages.require('pygrib')
    # from . import gribfile
except DistributionNotFound, exception:
    logging.warning('Python pygrib package is required for GRIB. ' 
        'No GRIB compatable modules are avaialble on this cerbere instance')
except ImportError:
    logging.exception('A GRIB mapper failed to load, and is unavailable:')

# from . import abstractmapper
# from . import urlseries

