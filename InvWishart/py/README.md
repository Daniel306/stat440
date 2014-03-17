Inverse Wishart Distribution in Python+Cython+Numpy+Scipy
=========================================================



Scrap Notes:
===========

Cython is a superset of python which allows you to--optionally, and only where needed--write C-ish code (real bare-metal for-loops, real machine ints and floats, etc). It evolved from pyrex which now has an unmaintained website.

'cython' the program is a compiler which takes .pyx files and writes them into .c files including all the necessary python-related #include statements
you then need to compile the .c file using `gcc` or `llvm` or whatever (making sure to -l and -I the correct libraries and #include directories) to a .so
the .so file is directly importable in Python, using python's built-in extension module abilities.


On my system, building under python2 creates "primes.so" and building under python3 gives "primes.cpython-33m.so"; this seems to be standard:
 except for PyQt4, sip, dbus, and the systemd bindings all python3 extension modules have this filename pattern, and all python2 do not.
