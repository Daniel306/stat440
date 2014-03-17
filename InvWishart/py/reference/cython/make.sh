# build primes python extension module, by hand
# --> this is just for reference; setuptools (somehow) can do this automagically
cython primes.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /usr/include/python3.3m/ -o primes.so primes.c
