#!/usr/bin/env python

from setuptools import setup


# Modified from http://stackoverflow.com/questions/2058802/
# how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package
def version():
    import os
    import re

    init = os.path.join('midtools', '__init__.py')
    with open(init) as fp:
        initData = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]+)['\"]",
                      initData, re.M)
    if match:
        return match.group(1)
    else:
        raise RuntimeError('Unable to find version string in %r.' % init)


setup(name='midtools',
      version=version(),
      packages=['midtools'],
      url='https://github.com/acorg/midtools',
      download_url='https://github.com/acorg/midtools',
      author='Terry Jones',
      author_email='tcj25@cam.ac.uk',
      keywords=['genetic sequences'],
      classifiers=[
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      description=('Command line scripts and Python classes for detecting '
                   'multiple infections from NGS reads.'),
      long_description='See https://github.com/acorg/midtools for details.',
      license='MIT',
      scripts=[
          'bin/base-frequencies.py',
          'bin/connected-components.py',
          'bin/consistency-basic.py',
          'bin/consistency-heatmap.py',
          'bin/coverage-and-significant-locations.py',
          'bin/create-mid-experiment-data.py',
          'bin/create-reads.py',
          'bin/multiple-significant-base-frequencies.py',
          'bin/mutate-reads.py',
          'bin/random-nt-sequence.py',
          'bin/significant-base-frequencies.py',
          'bin/simulation-significant-base-frequencies.sh',
      ],
      install_requires=[
          'dark-matter>=4.0.6',
          'plotly>=2.7.0',
          'flake8>=3.5.0',
      ])
