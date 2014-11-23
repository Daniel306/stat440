
from scipy import random
from numpy import *
from numpy.linalg import cholesky as chol, inv

Mu = array([1,3])
V = array([[148.84, 142.74], [142.74, 490.63]])
A = array([[3.38, -.77], [-.77, 2.55]])

def multigibbs(n, Mu, Sigma, j, thin=1, burnin=1):
	"""
	take n samples of
	  Y ~ N(Mu, Sigma)
	by gibbs-sampling the partitioned conditionals
	Y[:j] | Y[j:] and Y[j:] | Y[:j]
	Error checking is not done on arguments! in particular, make sure:
	- sigma is positive definite (and therefore symmetric)
	- mu and sigma have corresponding dimensions
	- 0 <= j < len(Mu)
	"""
	
	def g():
	
		d = len(Mu) #dimensionality...
		assert Mu.shape == (d,), "Mu must be a vector"
		assert Sigma.shape == (d,d), "Sigma must be a square matrix"
		assert (Sigma.T == Sigma).all(), "and symmetric"
	
		Mu_1, Mu_2 = Mu[:j], Mu[j:]
		Sigma_11, Sigma_12, Sigma_22 = Sigma[:j, :j], Sigma[:j, j:], Sigma[j:, j:]
	
		Y = array([0]*d) #it doesn't matter what we init Y to except that it is the right size, so init it to 0
	
		for i in range(n+burnin):
			for _ in range(thin): #skipe
				#TODO: precompute factors
				# also, you can probably get by faster by forcing Mu=0 during the gibbs sampling since then you avoid having to immediately subtract it off in the next tep, and then doing a single vectorized Y+=Mu at the end (or yield Y+Mu is close enough)
				Y[:j] = Mu_1 + \
					         dot(dot(Sigma_12, inv(Sigma_22)), Y[j:]-Mu_2) +\
					         dot(Sigma_11 - dot(dot(Sigma_12, inv(Sigma_22)), Sigma_12.T),
					             random.normal(size=j))
					             
				Y[j:] = Mu_2 + \
					         dot(dot(Sigma_12.T, inv(Sigma_11)), Y[:j]-Mu_1) +\
					         dot(Sigma_22 - dot(dot(Sigma_12.T, inv(Sigma_11)), Sigma_12),
					             random.normal(size=(d-j)))
			if i>=burnin:
				yield Y
				
				
	return array(list(g()))

def multinormal(n, Mu, Sigma):
	return random.multivariate_normal(Mu, Sigma, n)
