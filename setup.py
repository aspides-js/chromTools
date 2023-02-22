#!/usr/bin/env python

from setuptools import setup

exec(open("entirety/constants.py").read())

setup(name='entirety',
      version=COMPLETE_VERSION,
      package_dir = {'entirety' : 'entirety'},
      packages=['entirety'],
      scripts=['bin/ENTIRETY']
      )