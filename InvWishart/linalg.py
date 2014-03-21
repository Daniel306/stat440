"""
linalg.py

exploring the runtimes of various matrix operations by number crunching
"""


def matrix_multiply(A, B):
    "naive matrix multiply"
    "compute how long it would take to multiply matrix A and B"
    "don't actually pass a matrix; just A = (m,n) and B = (n, l)"
    m, n = A
    n2, l = B
    assert n == n2, "Incompatible matrix dimensions"
    
    S = {}  #our output matrix: a dictionary (just to be super inefficient)
    RUNTIME = 0  #TODO: break this into a dict counting muls, adds, and divides
    for i in range(m): #rows of the output matrix
        for j in range(l): #columns of the output matrix
            S[(i,j)] = 0
            for k in range(n):
               RUNTIME += 1
               # compute
    return RUNTIME*1.1 


# specializations on dense matrices:
# upper triangle*upper triangle
# lower triangle*lower triangle (can use the same alg as above; just need to transpose before and after)
# forward solve
# back solve (again, same alg)
# inversion (in general: requires backsolving
# inverting a triangle
# inve

# specializations on spare matrices:
# school book algorithm, but with skipping missing data

D = list(range(32))
naive_runtimes = [(d, matrix_multiply((d,d), (d,d))) for d in D]

import matplotlib.pyplot as plt
plt.plot(naive_runtimes, label="runtimes")
plt.plot(D, [d**2 for d in D], label="d^2")
plt.plot(D, [d**3 for d in D], label="d^3")
plt.legend()
plt.show()
