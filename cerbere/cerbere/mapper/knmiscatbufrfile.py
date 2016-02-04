# -*- coding: utf-8 -*-
"""
.. module::cerbere.mapper.knmiscatbufrfile

Mapper classs for KNMI L2 scatterometer BUFR format.

Requires the following packages:
  * https://github.com/pytroll/python-bufr

:copyright: Copyright 2013 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
.. codeauthor:: Jeff Piolle <jfpiolle@ifremer.fr>
"""

from cerbere.mapper import bufrfile


class KNMIScatBUFRFile(bufrfile.BUFRFile):
    """Mapper class to read KNMI L2 scatterometer BUFR files"""

    ATTRIBUTE_ENTRIES = {
        'satellite_identifier': 'platform',
        'cross_track_resolution': '',
        'along_track_resolution': '',
        'orbit_number': '',
        }
