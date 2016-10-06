import os
from setuptools import setup

# with open('README.rst') as readme_file:
#     README = readme_file.read()

# with open('HISTORY.rst') as history_file:
#     HISTORY = history_file.read()

VERSION = '0.1.0dev4'

LONG = ("Python 3 module implementing the Scanlon and Sahu (2008) procedure "
        "for partitioning measured water vapor and carbon dioxide fluxes into "
        "stomatal (transpiration, photosynthesis) and nonstomatal "
        "(evaporation, respiration) components.")

setup(
    name='fluxpart',
    version=VERSION,
    description='Partition water vapor and carbon dioxide fluxes',
    long_description=LONG,
    url='https://github.com/usda-ars-ussl/fluxpart',
    author='Todd Skaggs',
    author_email='todd.skaggs@ars.usda.gov',
    license='Public Domain CC0',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    packages=['fluxpart'],
    install_requires=['numpy', 'scipy', 'pywavelets'],
    zip_safe=False,
    tests_require=['pytest'],
    test_suite='tests',
)
