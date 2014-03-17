
from distutils.core import setup
from distutils.extension import Extension

setup(name="primes_demo",
     ext_modules = [Extension("primes", ["primes.pyx"])]
     )
