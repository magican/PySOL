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
    elif ext == '.zip':
        msg = """Files have to be extracted. \n Example:
        from zipfile import ZipFile
        import tempfile
        import shutil
        from cerbere.mapper.safegeotifffile import SAFEGeoTiffFile
        try:
            tmpFldr = tempfile.mkdtemp()  # create tmpdir
            with ZipFile(inFilename, "r") as zf:
                for name in zf.namelist():
                    localFilePath = zf.extract(name, tmpFldr)
                    #~ print localFilePath
            zf.close()
            sarmp = SAFEGeoTiffFile(url=localFilePath+filename)
        finally:
            try:
                shutil.rmtree(tmpFldr)  # delete directory
            except OSError as exc:
                if exc.errno != errno.ENOENT:  # ENOENT - no such file or directory
                    raise  # re-raise exception """
        raise Exception(msg)
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
