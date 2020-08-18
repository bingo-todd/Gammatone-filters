import os
from distutils.core import setup, Extension
os.environ['cc'] = "gcc"
setup(name='GTF',
      version='1.0',
      author='Todd',
      author_email='todd_stan@163.com',
      packages=['GTF'],
      package_data={'GTF': ['GTF.c']})
