import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name='meshcut',
    version='0.2.1',
    description='Utilities to slice 3D triangular meshes.',
    long_description=read('README.rst'),
    author='Julien Rebetez',
    author_email='julien@fhtagn.net',
    url='https://github.com/julienr/meshcut',
    download_url='https://github.com/julienr/meshcut/tarball/0.1',
    keywords=['mesh', 'slice', 'cross-section', '3D', 'triangular'],
    install_requires=[
        'numpy-stl',
        'scipy',
        'numpy',
    ],
    py_modules=['meshcut'],
)
