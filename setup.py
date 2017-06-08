#!/usr/bin/env python

from setuptools import setup, find_packages

config = {
    'name': 'stomasimulator',
    'description': 'Perform biomechanical simulations of stomata',
    'long_description': open('README.md').read(),
    'author': 'Hugh C. Woolfenden',
    'author_email': 'hugh.woolfenden@jic.ac.uk',
    'url': 'https://github.com/woolfeh/stomasimulator',
    'download_url': 'https://github.com/woolfeh/stomasimulator',
    'install_requires': ['pytest'],
    'packages': find_packages(),
    'version': '0.1.0',
    'scripts': ['bin/stomasimulator'],
    'license': open('LICENSE').read(),
    'include_package_data': True
}

setup(**config)
