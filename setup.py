#!/usr/bin/env python
from setuptools import setup

setup(name='Glutton',
      version='0.1',
      description='Transcriptome Scaffolding using Evolutionary Alignment',
      author='Alan Medlar',
      author_email='alan.j.medlar@helsinki.fi',
      url='https://upload.wikimedia.org/wikipedia/commons/b/b3/Gluttony.jpg',
      license='GNU Public License ver3 ( https://www.gnu.org/licenses/gpl-3.0.html )',
      long_description='Transcriptome Scaffolding using Evolutionary Alignment',
      platforms=['*nix'],
      packages=['glutton'],
      install_requires=['biopython>=1.6', 'sqlalchemy', 'pysam'],
      scripts=['scripts/glutton'],
     )
