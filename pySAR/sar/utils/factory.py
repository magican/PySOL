#!/usr/bin/env python
# coding=utf-8
"""
"""


import os


def sarmapper(filename):
    """Cerbere mapper factory."""
    ext = os.path.splitext(filename)[1].lower()
    if ext == '.tiff':
        from cerbere.mapper.safegeotifffile import SAFEGeoTiffFile
        sarmp = SAFEGeoTiffFile(url=filename)
    else:
        raise Exception('Unknown filename pattern.')
    return sarmp

def sarimage(filename):
    """
    """
    from sar.data.sarimage import SARImage
    sarmp = sarmapper(filename)
    sarim = SARImage(sarmp)
    return sarim
