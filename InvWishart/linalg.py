"""
linalg.py

exploring the runtimes of various matrix operations by number crunching.

Generic setup: two matrices of sizes (m,n) and (n,l).
Both are dense: all entries are (or at least, can't be ruled out) filled with nonzero values.

Basic matrix multiplication is at least O(m*l) (ie quadratic) because it must at least look at each index of the output matrix
 and the schoolbook way is O(m*n*l) (ie cubic) because inside each output index it does a dot product, which is O(n)
 and the best improvements we have on that use deep algebra cleverness which gives them absurdly huge hidden runtime constants.

  for special cases where the matrices are square (which is the case that can be extended to all other cases),
  I'll be loose and define d := m = n = l

# variant algorithms on dense matrices:
# upper triangle*upper triangle
# lower triangle*lower triangle (can use the same alg as above; just need to transpose before and after)
# forward solve (on a triangular matrix)
# backsolve (again, same alg as forwardsolve but with transposing before and after--which can be constant time, if you use numpy's array datastructure)
  - against a vector  (d^2: proof [...])
  - against a full matrix (same as as a vector, times a factor of d)
  - against a triangular matrix (of the same sidedness) (????)
  - against a triangular (of the opposite sidedness) (same as going against [...proof...])
# inversion (in general very expensive and also numerically unstable; the best algorithm is to LU decompose and then do a forwardsolve on the L and then a backsolve on the U
# inverting a triangle (
  - naive: backsolve against the identity matrix
  - more clever: only compute those points we know will e
# inver

# variants on spare matrices:
# school book algorithm, but with skipping missing data

We're only exploring dense algorithms here, because InvWishart depends on generating a Bartlett Factor which is dense in one triangle.


"""


def naive_dense_multiply(A, B):
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
                RUNTIME += 1; #for the multiply A[i,k] * B[k, j]
                RUNTIME += 1; # for the addition S[(i,j)] += {the above}
            # S[(i,j)] /= A[(i,j)]
            RUNTIME += 1 # for the division
    return RUNTIME 



D = list(range(32))
naive_runtimes = [(d, matrix_multiply((d,d), (d,d))) for d in D]

import matplotlib.pyplot as plt
plt.plot(naive_runtimes, label="naive dense multiply")
plt.plot(D, [d**2 for d in D], label="d^2")
plt.plot(D, [d**3 for d in D], label="d^3")
plt.legend()
plt.show()
