from setuptools import setup

with open('README.rst') as readme_file:
    README = readme_file.read()

with open('HISTORY.rst') as history_file:
    HISTORY = history_file.read()

setup(
    name='fluxpart',
    version='0.1.0dev1',
    description='Module for partitioning eddy covariance flux measurements.',
    long_description=README + '\n\n' + HISTORY,
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
