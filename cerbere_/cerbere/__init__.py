# encoding: utf-8
__all__ = ['datamodel', 'geo', 'mapper', 'plot', 'science', 'utils']

# Mapper constants, taken from abstractmapper to resolve cirular imports.
READ_ONLY = 'r'
WRITE_NEW = 'w'
READ_WRITE = 'r+'

# Status of mapper 
SAVED = 0
NOTSAVED = 1

# Module level default values
DEFAULT_TIME_UNITS = 'seconds since 1981-01-01 00:00:00'
CF_AUTHORITY = 'CF-1.8'

# from . import datamodel
# from . import geo
# from . import mapper
# from . import plot
# from . import science
# from . import utils

