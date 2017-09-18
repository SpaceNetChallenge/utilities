#!/usr/bin/env/python
"""Installation script
"""

import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


LONG_DESCRIPTION = """This repository has three python packages, geoTools and evalTools and labelTools. The geoTools packages is intended to assist in the preprocessing of [SpaceNet](https://spacenetchallenge.github.io/) satellite imagery data corpus hosted on [SpaceNet on AWS](https://aws.amazon.com/public-datasets/spacenet/) to a format that is consumable by machine learning algorithms. 
The evalTools package is used to evaluate the effectiveness of object detection algorithms using ground truth.
The labelTools package assists in transfering geoJson labels into common label schemes for machine learning frameworks
This is version 3.0 and has been updated with more capabilities to allow for computer vision applications using remote sensing data
"""

#if os.environ.get('READTHEDOCS', False) == 'True':
#    INSTALL_REQUIRES = []
#else:

INSTALL_REQUIRES = ['cython','geopandas', 'numpy', 'rtree', 'scipy', 'osmnx', 'centerline', 'affine', 'tqdm', 'rasterio>=1.0a9',
                    'pillow']

# get all data dirs in the datasets module
data_files = []

#for item in os.listdir("geopandas/datasets"):
#    if not item.startswith('__'):
#        if os.path.isdir(os.path.join("geopandas/datasets/", item)):
#            data_files.append(os.path.join("datasets", item, '*'))
#        elif item.endswith('.zip'):
#            data_files.append(os.path.join("datasets", item))


setup(name='spacenetutilities',
      version='3.0',
      description='Geographic pandas extensions',
      license='APACHE 2',
      author='SpaceNet Contributors',
      author_email='dlindenbaum@iqt.org',
      url='https://github.com/SpaceNetChallenge/utilities',
      long_description=LONG_DESCRIPTION,
      packages=['spacenetutilities',
                'spacenetutilities.labeltools',
                'spacenetutilities.scripts',
                ],
      package_data={'geopandas': data_files},
      install_requires=INSTALL_REQUIRES
      )
      #cmdclass=versioneer.get_cmdclass())