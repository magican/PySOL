#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cerbere.mapper.cryosat2ncfile
=============================

Mapper for CryoSat-2 Altimeter files from NOAA RADS 4.0

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""
from datetime import datetime

from .ncfile import NCFile


class Cryosat2NCFile(NCFile):
    """
    Mapper for CryoSat-2 Altimeter files from NOAA RADS 4.0.

    Overrides the NCFile mapper to take into account some specific
    attributes naming
    """

    def get_start_time(self):
        """Returns the minimum date of the file temporal coverage"""
        handler = self.get_handler()
        attrdate = handler.getncattr('first_meas_time')
        return datetime.strptime(attrdate, "%Y-%m-%d %H:%M:%S.%f")

    def get_end_time(self):
        """Returns the maximum date of the file temporal coverage"""
        handler = self.get_handler()
        attrdate = handler.getncattr('last_meas_time')
        return datetime.strptime(attrdate, "%Y-%m-%d %H:%M:%S.%f")
