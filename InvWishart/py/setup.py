
from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext as build_cython

# include_dirs = [numpy.get_include()]) <-- will need this

setup(name="primes_demo",
     packages = [''],
     cmdclass = {'build_ext': build_cython},
     ext_modules = [Extension("primes", ["primes.pyx"])]
     )
