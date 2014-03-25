Proofs
=======

We require derivations (or maybe refutations) of these facts:

* Bartlett's Theorem: if A[i,i] ~ sqrt(X^2_{df-{i-1}}) && A[i,j] ~ N(0,1) and A[j,i] = 0 (i>j) then AA' ~ W(I, df)
  * The Wishart linearity lemma: if M ~ W(V,df) then CMC' ~ W(CVC', df) (but only when squared like this and when C is full fank)
  * the equivalent for inverse wisharts:
    * if M ~ W^-1(V, df), then CMC' ~ W^-1(???, df)
* Facts about symmetric matrices?
* Facts about positive definite matrices
* Fact about the Cholesky Decomposition and Upper/Lower Triangles
  * Uniqueness(??)
  * That it can be thought of in two ways: lower-triangles or upper triangles
  * Algorithmic speed (d^3/3, half less than generic LU )
  * Relationship between the LU (cholesky) and UL (<-- this isn't anywhere???) decompositions (and the relationship to inversion)
  * Invertibility of cholesky factors: rank(LU) = min(rank(L), rank(U)) ={in the cholesky decomp L = U'} = min(rank(L), rank(L')) = rank(L)
  * if chol(M) is the right term of the decomposition then t(chol(M)) is the left
  * chol(inv(M)) = inv(chol(M))
  * the UL decomposition: U = solve(chol(solve(M))) ==> UU' = M (upper times a lower!)
* Lower * Lower = Lower and Upper*Upper = Upper but Lower*Upper = anything ((so lower*t(lower) = anything))
* inv(Lower) = Lower (corollary: inv(Upper) = Upper, since upper = lower' and inverse swaps with transposing: inv(Lower)' = Lower' = Upper = inv(Lower') = inv(Upper)       
* solve(backsolve(M, solve(t(M)))) = M'M

* In HW1 we showed that the mean of the Wishart distribution is the covariance matrix of the associated Normals (scaled by a constant: df).
  In terms of types, that means the Wishart distribution has type "covariance matrix". And that should mean that the Inverse-Wishart has type "https://en.wikipedia.org/wiki/Precision_matrix".  All of the references I've seen elsewhere seem to treat the output of the Inverse-Wishart as another covariance matrix. Why is this?
  What is up with this?
