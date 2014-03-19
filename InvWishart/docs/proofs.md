Proofs
=======

We require derivations (or maybe refutations) of these facts:

* Bartlett's Theorem: if A[i,i] ~ sqrt(X^2_{df-{i-1}}) && A[i,j] ~ N(0,1) and A[j,i] = 0 (i>j) then AA' ~ W(I, df)
* The Wishart linearity lemma: if M ~ W(V,df) then CWC' ~ W(CVC', df) (but only when squared like this)
* Facts about symmetric matrices
* Facts about positive definite matrices
* Fact about the Cholesky Decomposition and Upper/Lower Triangles
  * Uniqueness(??)
  * That it can be thought of in two ways: lower-triangles or upper triangles
  * Algorithmic speed (d^3/3, half less than generic LU )
  * Relationship between the LU (cholesky) and UL (<-- this isn't anywhere???) decompositions (and the relationship to inversion)
  * Invertibility of cholesky factors: rank(LU) = min(rank(L), rank(U)) = rank(L) because rank(L) = rank(L')
  * if chol(M) is the right term of the decomposition then t(chol(M)) is the left
  *  chol(inv(M)) = t(inv(chol(M)))     (corollary: t(chol(inv(M))) = inv(chol(M)), so you can swap inversions inside the chol for transposes and inverses outside it or vice-versa)
     and inv(chol(inv(M))) = chol(M)'
  * the UL decomposition: U = solve(chol(solve(M))) ==> UU' = M (upper times a lower!)
  * chol(M
* Lower * Lower = Lower and Upper*Upper = Upper but Lower*Upper = anything ((so lower*t(lower) = anything))
* inv(Lower) = Lower (corollary: inv(Upper) = Upper, since upper = lower' and inverse swaps with transposing: inv(Lower)' = Lower' = Upper = inv(Lower') = inv(Upper)       
