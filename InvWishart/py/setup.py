
from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext as build_cython

setup(name="primes_demo",
     cmdclass = {'build_ext': build_cython},
     ext_modules = [Extension("primes", ["primes.pyx"])]
     )
