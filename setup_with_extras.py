#!/usr/bin/env python
from setuptools import setup

setup(name='Glutton',
      version='0.1',
      description='Transcriptome scaffolding and postprocessing for comparative analysis using evolutionary alignment',
      author='Alan Medlar',
      author_email='alan.j.medlar@helsinki.fi',
      url='http://wasabiapp.org/software/glutton',
      license='GNU Public License ver3 ( https://www.gnu.org/licenses/gpl-3.0.html )',
      long_description='Transcriptome scaffolding and postprocessing for comparative analysis using evolutionary alignment',
      platforms=['*nix'],
      packages=['glutton'],
      install_requires=['biopython>=1.6', 'sqlalchemy', 'mysql-python', 'pysam'],
      scripts=['scripts/glutton', 'binaries/prank', 'binaries/pagan', 'binaries/exonerate', 'binaries/blastx', 'binaries/bppphysamp'],
     )

