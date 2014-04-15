def primes(int kmax):
    cdef int n, k, i
    cdef int p[1000]
    result = []
    if kmax > 1000:
        kmax = 1000
    k = 0
    n = 2
    while k < kmax:
        i = 0
        while i < k and n % p[i] != 0:
            i = i + 1
        if i == k:
            p[k] = n
            k = k + 1
            result.append(n)
        n = n + 1
    return result



print(__name__)
if __name__ == '__main__':
    # doesn't run. this module needs to be compiled into an extension module (cython .pyx -> gcc .c -> pythin import .so)
    # and as an extension module, it can never be run directly
    print(primes(55))
