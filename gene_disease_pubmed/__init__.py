"""PubMed abstract query post processor.


Copyright (c) 2014 Jeroen F.J. Laros

Licensed under the MIT license, see the LICENSE file.
"""
import argparse
import os


__version_info__ = ('0', '0', '1')

__version__ = '.'.join(__version_info__)
__author__ = 'Jeroen F.J. Laros'
__contact__ = 'J.F.J.Laros@lumc.nl'
__homepage__ = 'https://git.lumc.nl/j.f.j.laros/gene-disease-pubmed/tree/master'

usage = __doc__.split('\n\n\n')


class ProtectedFileType(argparse.FileType):
    def __call__(self, string):
        if 'w' in self._mode and os.path.exists(string):
            raise IOError('failed to create "{}": file exists.'.format(string))
        return super(ProtectedFileType, self).__call__(string)


def doc_split(func):
    return func.__doc__.split('\n\n')[0]


def version(name):
    return '{} version {}\n\nAuthor   : {} <{}>\nHomepage : {}'.format(
        name, __version__, __author__, __contact__, __homepage__)
