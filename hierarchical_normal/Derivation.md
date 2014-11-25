Derivation

1) define \sigma = chol(\Sigma), noting that \sigma is a) lower triangular b) you can write Y = U + \sigma Z to get Y ~ N(U, \Sigma)
2) partition Y, \sigma and Z
3) work out what the partitioned terms 


sigma:
(note: x^2 means tcrossprod(x), i.e. in numpy notation dot(x,x.T) and is the inverse operation to chol() which computes the LU decomposition--at least, inverse on the space of things which can be chol'd)
namely:
Sigma_11 = sigma_11^2
Sigma_12 = sigma_11 sigma_21^T 
Sigma_22 = sigma_21^2 + sigma_22^2
Sigma_21 = Sigma_21^T
 finally
0        = sigma_12, since sigma is chol(Sigma) which makes it *lower triangular*. This part is key! Without it we cannot do the conditioning--not without dropping down to integrating PDFs, which sucks.


4) from that + some algebra, derive these two distributions:
Y_1 ~ N(U_1, \sigma_11 * \sigma_11^T)
Y_2 | Y_1 ~ N(sigma_21 inv(sigma_11)(Y_1 - U_1) + U_2, \sigma_22 * \sigma_22^2)

then, toss out \sigma in favour of \Sigma as much as possible:
sigma_11^2 = Sigma_11
sigma_21 inv(sigma_11) = Sigma_12^T inv(sigma_11)^T inv(sigma_11)
                       = Sigma_12^T inv(Sigma_11)
\sigma_21^2 has a very useful form:
 = (Sigma_12^T inv(sigma_11)^T) (inv(sigma_11) Sigma_12)
 = Sigma_21 inv(sigma_11 sigma_11^T) Sigma_12
 = Sigma_21 inv(Sigma_11) Sigma_12
\sigma_22^2 = \Sigma_22 - \sigma_21^2
            = Sigma_22 - Sigma_21 inv(Sigma_11) Sigma_12

So    Y_1 ~ N(U_1, Sigma_11)
Y_2 | Y_1 ~ N(Sigma_21 inv(Sigma_11)(Y_1 - U_1) + U_2, Sigma_22 - Sigma_21 inv(Sigma_11) Sigma_12)

the important result is:
Y1|Y2 ~ N(Sigma_12 inv(Sigma_22)(Y2 - U_2) + U_1,
          Sigma_11 - Sigma_12 inv(Sigma_22) Sigma_21)
and by symmetry--just reorder all the rows ( Y = [Y1; Y2], and the corresponding reordering of Sigma as well)--you can get 
Y2|Y1 ~ N(Sigma_21 inv(Sigma_11)(Y1 - U_1) + U_2,
          Sigma_22 - Sigma_21 inv(Sigma_11) Sigma_12)

I don't know a good name for this. I've been calling it the conditional multinormal formula.
 a description of it, without proof, but with applications, is at <>

---

then, to actually apply to the hierarchical normal, define:

we desire a pair of linked multinormal distributions such that
    U ~ N(XB, A)
Y | U ~ N(U, V)


((( note that this is equivalent to doing
which a = chol(A), v = chol(V)

U = XB + aZ1
Y = U + vZ2
 TODO: work out how to have these formulas fall out of the structure above)

like, show that 
Y2 = U_2 + sigma_21*z1 + sigma_22*z2
   = Y1 + sigma_22*z2
)))


(the idea being that you A and B (and maybe V, tho you might just set that to I or something simple--it's sort of irrelevant) are model parameters, each U is an intermediate "random effect":
 Y = XB + random_effect + error
 you may ask "why can't the random effect be counted as error?"
  the answer is that the random effects are uncorrelated and possibly have different variances, so homoscedasticity fails
 e.g. maybe each sample comes from 
 (more usefully, model it as many samples associated with a single random effect, e.g. samples = students, "effect" = school)
we can do this with the above by setting
U = Y_1 
Y = Y_2

XB = U_1
A = Sigma_11
U = Sigma_21 inv(Sigma_11)(U - U_1) + U_2
V = Sigma_22 - Sigma_21 inv(Sigma_11) Sigma_12

we need to work out what U_2, Sigma_21, and Sigma_22 are.
 If we know Sigma_21 then we know Sigma_12 and therefore can solve the last equation for Sigma_22,
 and similarly can solve the second last for U_2
 so the sticking point is Sigma_21

so, get Sigma_21 directly by its definition: it is the matrix whose entries are the covariance of each individual entry
Sigma_12 = cov(U, Y)
         = E[(U - E[U])(Y - E[Y])^T]
         = E[(XB + aZ1 - E[XB + aZ1]) (XB + aZ1 + vZ2 - E[XB + aZ1 + vZ2])^T]
         = E[(aZ1) (aZ1 + vZ2)^T]
         = E[(aZ1Z1^Ta^T + aZ1Z2^Tv^T]
         = aE[Z1 Z1^T]a^T + aE[Z1 Z2^T]v^T  #<-- first term = var(Z1) since E[Z1] = 0. second term = cov(Z1,Z2) = E[Z1Z2] - E[Z1]E[Z2] = E[Z1Z2] - 0*0. and we know Z1 and Z2 are independent so this is 0
         = aIa^T  + a0b^T
         = aa^T
         = A
         = Sigma_11

so Sigma_22 = V + A A^-1 A = V + A
so immediately we know the full multivariate normal covariance is Sigma = [A A; A A+V]

and
U = A inv(A)(U - XB) + U_2
U = U - XB + U_2
XB = U_2
 --> U_1 = U_2 ((which makes a lot of intuitive sense: the mean of the second random variable has to be the same as the first, because it uses the first as its (conditional) mean))

((
 if this result seems strange since Sigma_21 should somehow be fundamentally of different character than Sigma_11, remember that in this particular case, the partitioning is done exactly down the middle.
  [also i have experimental support that this is the right formula, see multigibbs.py, so take that]
))

and the *conditionals*, by combining this result and the result above

Y | U ~ N(Sigma_21 inv(Sigma_11)(Y1 - U_1) + U_2,
          Sigma_22 - Sigma_21 inv(Sigma_11) Sigma_12)
       =N(A inv(A)(U - U_1) + U_2,
          (A+V) - A inv(A) A)
       =N(U - U_1 + U_1, A + V - A)
       =N(U, V)
        as desired and asserted
   and then for the prize as far as gibbs sampling is concerned: the *reverse* conditional: 
       
U | Y ~ N(Sigma_12 inv(Sigma_22)(Y2 - U_2) + U_1,
          Sigma_11 - Sigma_12 inv(Sigma_22) Sigma_21)
       =N(A inv(A+V)(Y - U_2) + U_1,
          A - A inv(A+V) A)
       =N(A inv(A+V)(Y - XB) + XB,   A(I-inv(A+V)A))
       =N(A inv(A+V)(Y - XB) + XB, A(I-inv(A+V)*A))


 
------------


Curious:
 the notes (between typos) give a different formula for the reverse:
    U | Y ~ N(V inv(V+A)(XB - Y) + Y, (I-V*inv(V+A))*V) 
 which is so similar as to make you think something is up
 but also just slightly different
 
 experimentally, this formula seems to work, so take that for what it's worth