# -*- coding: utf-8 -*-
"""
cerbere.datamodel.variable
==========================

Classes for defining a geophysical quantity

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

__all__ = ['Variable']


class Variable(object):
    """Description of the phenomenon (or measured quantity) related to a
    measure (temperature, salinity,...)
    """
    def __init__(self, shortname, description=None,
                 authority=None, standardname=None):
        """
        Args:
            shortname (string): the label of the variable (don't use any white
                space). This corresponds to the variable name in a netcdf file

            description (string) : full name of the phenomenon. This corresponds
                to a long_name in attribute in a netCDF file.

            authority (string): naming authority referencing the provided
                standard name

            standardname (string): standard label for a phenomenon, with
                respect to the convention stated in `authority` argument.
                This corresponds to a standard_name attribute in a CF compliant
                NetCDF file.
        """
        object.__init__(self)
        self.shortname = shortname
        self.description = description
        self.authority = authority
        self.standardname = standardname
