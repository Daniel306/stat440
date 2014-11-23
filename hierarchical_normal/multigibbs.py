
from scipy import random
from numpy import *
from numpy.linalg import cholesky as chol, inv

import matplotlib.pyplot as plt

Mu = array([1,3])
V = array([[148.84, 142.74], [142.74, 490.63]])
A = array([[3.38, -.77], [-.77, 2.55]])

def multigibbs(n, Mu, Sigma, j, thin=10, burnin=5000):
	"""
	take n samples of
	  Y ~ N(Mu, Sigma)
	by gibbs-sampling the partitioned conditionals
	  Y[:j] | Y[j:] and
	  Y[j:] | Y[:j]
	
	This is just a prototype! 'j' should not be exposed as part of the API.
	"""
	
	# To gibbs-sample you need to know a chain of conditional distributions
	# In this case, we need two: Y1 | Y2 and Y2 | Y1
	# for which we have the conditional multinormal formula
	# e.g. https://onlinecourses.science.psu.edu/stat505/node/43
	# Y1 | Y2 ~ N(mu1 + cov(Y1,Y2).inv(var(Y2)).(Y2-mu2),  var(Y1) - cov(Y1,Y2).inv(var(Y2)).cov(Y2,Y1))
	# and vice versa for Y2 | Y1
	
	def g():
		
		d = len(Mu) #dimensionality...
		assert Mu.shape == (d,), "Mu must be a vector"
		assert Sigma.shape == (d,d), "Sigma must be a square matrix"
		assert (Sigma.T == Sigma).all(), "and symmetric"
		assert 0<=j<d
		
		Mu_1, Mu_2 = Mu[:j], Mu[j:]
		Sigma_11, Sigma_12, Sigma_22 = Sigma[:j, :j], Sigma[:j, j:], Sigma[j:, j:]
	
		Y = array([0.0]*d) #it doesn't matter what we init Y to except that it is the right size, so init it to 0
	
		for i in range(n+burnin):
			for _ in range(thin): #skipe
				#TODO: precompute factors
				# also, you can probably get by faster by forcing Mu=0 during the gibbs sampling since then you avoid having to immediately subtract it off in the next tep, and then doing a single vectorized Y+=Mu at the end (or yield Y+Mu is close enough)
				
				if j > 0:
					Y[:j] = Mu_1 + \
							     dot(dot(Sigma_12, inv(Sigma_22)), Y[j:]-Mu_2) +\
							     dot(chol(Sigma_11 - dot(dot(Sigma_12, inv(Sigma_22)), Sigma_12.T)),
							         random.normal(size=j))
				if j < d:
					Y[j:] = Mu_2 + \
							     dot(dot(Sigma_12.T, inv(Sigma_11)), Y[:j]-Mu_1) +\
							     dot(chol(Sigma_22 - dot(dot(Sigma_12.T, inv(Sigma_11)), Sigma_12)),
							         random.normal(size=(d-j)))
			if i>=burnin:
				yield array(Y)
				
				# reasons this is wrong:
				# - the sigmas should be chol()'d first.
				#   in the univariate normal, a*z makes a variable with a**2 variance
				#   the same applies when matrices happen
				# - the sigmas
				
				# edge cases to check:
				# j = 0 or j = d
				# that zero covariance is the same
				# that this gives the same results as m
				
				# TODO:
				# - to use this to sample the hierarchical normal,
				#   i need to work out the covariance of u and y
				#   ((lysy gave us the variances of each: V and A, so we know Sigma_11 and Sigma_22, but not Sigma_12))
				
	return array(list(g()))

def multinormal(n, Mu, Sigma):
	return random.multivariate_normal(Mu, Sigma, n)


def test_compare2d(D1, D2, i, j):

	plt.scatter(D1[:,i], D1[:,j])
	plt.title("Ground truth sample")
	plt.figure()
	plt.scatter(D2[:,i], D2[:,j])
	plt.title("Gibbs sample")
	plt.show()

def test_multigibbs():
	d = 18
	Mu = random.uniform(-9, 9, size=d)
	A = random.uniform(-33, 33, size=(d,d)); A = dot(A,A.T) #crap way of making a positive semidef matrix
	
	n = 99999
	Sd = multinormal(n, Mu, A)
	j = d//3 + 1
	Sg = multigibbs(n, Mu, A, j=j)
	
	try:
		# these should be approximately the same
		print("Expected mean:")
		print(Mu)
		print("Expected variance:")
		print(A)
		
		
		print("built-in multinormal sampler")
		print("mean (largest absolute deviation):")
		print(abs(Mu - Sd.mean(axis=0)).max())
		print("variance (largest absolute deviation):")
		print(abs(A - cov(Sd.T)).max())
	
		print("gibbs sampler at j=", j)
		print("mean (largest absolute deviation):")
		print(abs(Mu - Sg.mean(axis=0)).max())
		print("variance (largest absolute deviation):")
		print(abs(A - cov(Sg.T)).max())
	
		test_compare2d(Sd, Sg, 0, 1)
		#import IPython; IPython.embed() #DEBUG
	finally:
		return Sd, Sg #so that python -i can do fun things
	
if __name__ == '__main__':
	Sd, Sg = test_multigibbs()