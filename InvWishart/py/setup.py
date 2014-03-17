
from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext as build_cython

build_cython.inplace = True  #make the output come out somewhere obvious
                             #XXX this is incompatible with 'clean's default behaviour; I get 'primes.c' left uncleaned

# include_dirs = [numpy.get_include()]) <-- will need this eventually

setup(name="primes_demo",
     packages = [''], # <-- will need this eventually
     cmdclass = {'build_ext': build_cython},
     ext_modules = [Extension("primes", ["primes.pyx"])]
     )
