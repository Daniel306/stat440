
from scipy import random
from numpy import *
from numpy.linalg import cholesky as chol, inv

import matplotlib.pyplot as plt


def direct_hierarchical_normal(n, V, Mu, A):
	"""
	sample
		(y1,y2) such that
		y1 ~ N(Mu, A)
		y2|y1 ~ N(y1, V)
	directly from the definition
	(for comparison)
	"""
	a = chol(A)
	v = chol(V)
	d = len(Mu)
	
	def model():
		z1, z2 = random.normal(size=d), random.normal(size=d)
		y1 = Mu + dot(a, z1)
		y2 = y1 + dot(v, z2)
		return y1, y2
	
	return array([model() for i in range(n)])

def multinormal_hierarchical_normal(n, V, Mu, A):
	"""
	sample the hierarchical normal according to what I derived 
	"""
	d = len(Mu)
	r = multinormal(n, r_[Mu, Mu], r_[c_[A, A], c_[A, A+V]])
	r.shape = (n, d, d)
	return r

def gibbs_hierarchical_normal(n, V, Mu, A):
	"""
	sample the hierarchical normal according to what I derived...
	 but by falling back on the generic multinormal gibbs sampler
	"""
	d = len(Mu)
	r = multigibbs(n, r_[Mu, Mu], r_[c_[A, A], c_[A, A+V]], j=d)
	r.shape = (n, d, d)
	return r

def gibbs_hierarchical_normal_unrolled(n, V, Mu, A, thin=10, burnin=5000):
	"""
	sample the hierarchical normal according to what I derived *with unrolled formulas*
	((i suspect it is not actually any faster to do it this way))
	"""
	def g():
		d = len(Mu)
		assert Mu.shape == (d,), "Mu must be a vector"
		assert A.shape == (d,d), "A must be a square matrix"
		assert (A.T == A).all(), "and symmetric"
		assert V.shape == (d,d), "V must be a square matrix"
		assert (V.T == V).all(), "and symmetric"
		
		a = chol(A)
		v = chol(V)
		
		B = dot(A, inv(A+V))
		_a2 = A - dot(B, A)
		_a2 = chol(_a2)
		
		Y, U = array([0.0]*d), array([0.0]*d)
		
		for i in range(n+burnin):
			for _ in range(thin): #skipe
				# sample Y | U ~ N(U, V)
				Y = U + dot(v, random.normal(size=d))
				
				# sample U | Y ~ N(A(A+V)^-1)(Y-Mu) + Mu,
				#                  A - A(A+V)^-1A)
				U = dot(B, (Y-Mu)) + Mu +\
					 + dot(_a2, random.normal(size=d))
				
			if i>=burnin:
				yield [U, Y]
	
	return array(list(g()))
	
def gibbs_hierarchical_normal_wrong(n, V, Mu, A, thin=10, burnin=5000):
	"""
	sample the hierarchical normal according to what the assignment notes say
	 I *suspect* that the notes have A and V, and Mu and y2, swapped
	  (and also have a , typo'd where a + should be)
  	
  	this *should* come out with the wrong numbers
  	though maybe there's a hilarious third derivation that i missed
  	...huh. 
	"""
	def g():
		d = len(Mu)
		assert Mu.shape == (d,), "Mu must be a vector"
		assert A.shape == (d,d), "A must be a square matrix"
		assert (A.T == A).all(), "and symmetric"
		assert V.shape == (d,d), "V must be a square matrix"
		assert (V.T == V).all(), "and symmetric"
		
		a = chol(A)
		v = chol(V)
		
		B = dot(V, inv(V + A))
		_a2 = V - dot(B, V)
		_a2 = chol(_a2)  
		
		Y, U = array([0.0]*d), array([0.0]*d)
		
		for i in range(n+burnin):
			for _ in range(thin): #skipe
				# sample Y | U ~ N(U, V)
				Y = U + dot(v, random.normal(size=d))
				
				# sample U | Y ~ N(A(A+V)^-1)(Y-Mu) + Mu,
				#                  A - A(A+V)^-1A)
				U =  dot(B, (Mu-Y)) + Y +\
					 + dot(_a2, random.normal(size=d))
				
			if i>=burnin:
				yield [U, Y]
	
	return array(list(g()))
	

#----------------------------------------------------

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

def test_hierarchical_normal():	
	Mu = array([0,0]) #xi = 1, B=(0,0), which is a pretty lame test case to be honest
	V = array([[148.84, 142.74], [142.74, 490.63]])
	A = array([[3.38, -.77], [-.77, 2.55]])
	
	d, d = V.shape
	
	n = 2323
	
	# each sample is a 2x2 matrix:

	print("expected covariances:")
	print(A)
	print(A)
	print(A)
	print(A + V)
	print()
	
	for sampler in [direct_hierarchical_normal,
					multinormal_hierarchical_normal,
					gibbs_hierarchical_normal,
					gibbs_hierarchical_normal_unrolled,
					gibbs_hierarchical_normal_wrong,
					]:
		S = sampler(n, V, Mu, A)
		original_shape = S.shape
		
		# flatten so that we can use
		S.shape = (n, 2*d)	
		# 'cov', which gets confused on non-matrix input 
		C = cov(S.T)
		S.shape = original_shape
	
		assert C.shape == (2*d, 2*d)
		C.shape = (2, d, 2, d)
		
		print(sampler.__name__)
		for i in range(2):
			for j in range(2):
				print(C[i, :, j, :])
		print()
		
	import IPython; IPython.embed() #DEBUG
	
if __name__ == '__main__':
	#Sd, Sg = test_multigibbs()
	test_hierarchical_normal()