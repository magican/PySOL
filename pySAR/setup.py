# encoding: utf-8

"""
sar
===

Utilities to handle SAR data.
"""
from setuptools import setup

setup(
    zip_safe=False,
    name='sar',
    version='0.0.1',
    author='Gilles Guitton <gilles.guitton@oceandatalab.com>',
    author_email='gilles.guitton@oceandatalab.com',
    packages=[ 'sar'
             , 'sar.data'
             , 'sar.external'
             , 'sar.processing'
             , 'sar.render'
             , 'sar.transform'
             , 'sar.utils'
             , 'sar.utils.datapath'
             ],
    scripts=[
        'bin/s1-path',
        'bin/sar_crossspectra_png.py',
        'bin/sar_dn_png.py',
        'bin/sar_roughness_png.py',
        'bin/sar-xspec'
    ],
    license='LICENSE.txt',
    description='Utilities to handle SAR data.',
    long_description=open('README.txt').read(),
    install_requires=[ 'numpy'
                     , 'scipy'
                     , 'shapely'
                     , 'matplotlib'
                     , 'cerbere'
    ],
    package_data={ 'sar': [ 'share/logos/*.*'
                          , 'render/palette/*.*']
    },
)
