import os
import sys
import platform
from os import sep as dirsep
from os.path import isfile, join

from setuptools import setup, Extension
from Cython.Build import cythonize

if sys.version_info[:2] < (2, 6):
    sys.stderr.write('Python 2.5 and older is not supported\n')
    sys.exit()

if os.name == 'java':
    sys.stderr.write('JavaOS is not supported\n')
    sys.exit()

try:
    import numpy
except ImportError:
    sys.stderr.write('numpy is not installed, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

try:
    import prody
except ImportError:
    sys.stderr.write('prody is not installed, you can find it at: '
                     'http://http://prody.csb.pitt.edu\n')
    sys.exit()

if [int(dgt) for dgt in numpy.__version__.split('.')[:2]] < [1, 4]:
    sys.stderr.write('numpy v1.4 or later is required, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

__version__ = ''

with open('README.md') as inp:
    long_description = inp.read()

extensions = [
    Extension('dfsutils.libdfs', ['dfsutils/libdfs.pyx'],
    		include_dirs=[numpy.get_include()]) 
    ]

packages = ['dfsutils']

scripts = ['scripts/dfs']



setup(
    name='DFS',
    version=__version__,
    author='Matteo Tiberti',
    author_email='matteo.tiberti@gmail.com',
    description='A Python Package for the prediction of compensatory mutations in proteins',
    long_description=long_description,
    url='http://www.google.com',
    packages=packages,
    ext_modules=cythonize(extensions),
    license='GPL v3',
    keywords=('protein, dynamics, elastic network model, '
              'Gaussian network model, anisotropic network model, '
              'Protein Data Bank, PDB, ANM'),
    classifiers=[
                 'Development Status :: 5 - Production/Stable',
                 'Intended Audience :: Education',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: MIT License',
                 'Operating System :: MacOS',
                 'Operating System :: POSIX',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 2',
                 'Programming Language :: Python :: 3',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Chemistry',
                ],
    scripts=scripts
    #requires=['NumPy (>=1.7), ProDy'],
    #provides=['DFS ({0:s})'.format(__version__)]
)
