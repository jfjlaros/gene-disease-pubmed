import os
import sys
from setuptools import setup

requires = ['biopython', 'restkit']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

# This is quite the hack, but we don't want to import our package from here
# since that's recipe for disaster (it might have some uninstalled
# dependencies, or we might import another already installed version).
distmeta = {}
for line in open(os.path.join('gene_disease_pubmed', '__init__.py')):
    try:
        field, value = (x.strip() for x in line.split('='))
    except ValueError:
        continue
    if field == '__version_info__':
        value = value.strip('[]()')
        value = '.'.join(x.strip(' \'"') for x in value.split(','))
    else:
        value = value.strip('\'"')
    distmeta[field] = value

#try:
#    with open('README.md') as readme:
#        long_description = readme.read()
#except IOError:
long_description = 'See ' + distmeta['__homepage__']

setup(
    name='gene-disease-pubmed',
    version=distmeta['__version_info__'],
    description='PubMed abstract query post processor.',
    long_description=long_description,
    author=distmeta['__author__'],
    author_email=distmeta['__contact__'],
    url=distmeta['__homepage__'],
    license='MIT License',
    platforms=['any'],
    packages=['gene_disease_pubmed'],
    install_requires=requires,
    entry_points={
        'console_scripts': [
            'gdp = gene_disease_pubmed.gene_disease_pubmed:main'
        ]
    }
)
