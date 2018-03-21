"""
py2app/py2exe build script for ALPHA.

Will automatically ensure that all build prerequisites are available 
via ez_setup.

Usage (Mac OS X):
    python setup.py py2app

Usage (Windows):
    python setup.py py2exe
"""

from ez_setup import use_setuptools
use_setuptools()

import sys
from setuptools import setup

APP = ['source/main.py']
DATA_FILES = [
	'imgs/alphaLogo.png',
	'imgs/LStatisticTree.png',
	'imgs/warning.png',
	'imgs/tree.png',
]
PACKAGES = [
	'sip', 
	'PyQt4._qt', 
	're', 
	'os', 
	'itertools', 
	'ete3', 
	'copy', 
	'subprocess', 
	'collections', 
	'scipy', 
	'natsort', 
	'functools', 
	'webbrowser', 
	'os', 
	'matplotlib', 
	'cStringIO', 
	'numpy', 
	'math', 
	'random', 
	'shutil', 
	'dendropy', 
	'PIL', 
	'reportlab',
    'reportlab.lib',
    'reportlab.lib.colors',
    'reportlab.platypus',
    'biopython'
 ]
OPTIONS = {
	'iconfile': 'alphaLogo.icns', 
	'argv_emulation': True, 
	'includes': PACKAGES
}
# if sys.platform == 'darwin':
	# pass
# elif sys.platform = 'win32'
setup(
	name='ALPHA',
	version='1.0',
	description='Automated Local Phylogenomic Analyses',
	url='https://github.com/chilleo/ALPHA',
	author='Chabrielle Allen, Travis Benedict, Ryan Elworth, Peter Dulworth',
	author_email='chileo@gmail.com',
	license='MIT',
    app=APP,
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    install_requires=['pillow'],
    setup_requires=['py2app'],
)


#  mainscript = 'MyApplication.py'

#  if sys.platform == 'darwin':
#      extra_options = dict(
#          setup_requires=['py2app'],
#          app=[mainscript],
#          # Cross-platform applications generally expect sys.argv to
#          # be used for opening files.
#          options=dict(py2app=dict(argv_emulation=True)),
#      )
#  elif sys.platform == 'win32':
#      extra_options = dict(
#          setup_requires=['py2exe'],
#          app=[mainscript],
#      )
# else:
#      extra_options = dict(
#          # Normally unix-like platforms will use "setup.py install"
#          # and install the main script as such
#          scripts=[mainscript],
#      )

# setup(
#     name="MyApplication",
#     **extra_options
# )
