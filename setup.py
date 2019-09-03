from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name = 'ProtCHOIR',
  packages = ['ProtCHOIR'],
  version = '0.1.2',
  description = '[NOT RELEASED] A modeller-based pipeline to generate homo-oligomers.',
  long_description = long_description,
  license='CC-BY-4.0',
  author = 'Pedro Torres, Sony Malhotra, Tom Blundell',
  author_email = 'phm34@cam.ac.uk.com',
  url = 'https://github.com/monteirotorres/ProtCHOIR',
  download_url = 'https://github.com/monteirotorres/ProtCHOIR/archive/0.100.tar.gz',
  keywords = ['PDB',
              'structure',
              'protein',
              'oligomerization',
              'oligomers',
              'multimer'],
  classifiers = [],
  py_modules=[
      'ProtCHOIR.ProtCHOIR',
      'ProtCHOIR.Toolbox',
      'ProtCHOIR.Initialise'],
  entry_points={
      'console_scripts': [
          'ProtCHOIR = ProtCHOIR.__main__.py:main']},
  install_requires=[
      'markdown',
      'biopandas',
      'pandas',
      'progressbar2',
      'pathlib2',
      'biopython',
      'parasail',
      'matplotlib',
      'modeller',
      'pickle',
      'textwrap',
      'gzip',
      'numpy',
      'jinja2',
      'networkx'])
