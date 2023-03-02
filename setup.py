#!/usr/bin/env python

from setuptools import setup

exec(open("entirety/constants.py").read())

setup(name='entirety',
      version=COMPLETE_VERSION,
      packages=['entirety'],
      scripts=['bin/ENTIRETY'],
      package_data={'entirety': ['chromsize/*.txt']}
      )