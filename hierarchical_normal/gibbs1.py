"""
univariate prototyping of flipping conditioning
 in order to be able to do the multivariate gibbs sampler

"""

from numpy.random import normal

def gibbs1(burnin=10000, n=1000000, thin=50):
    """
    sample 
      u ~ N(3, 16)
    y|u ~ N(u, 25)
    this is a random effects model
    
    via gibbs sampling
    """
    
    u, y = 0,0
    
    for i in range(burnin+n):
        for _ in range(thin):
            y = u + 5*normal()
            u = y - 5*normal()
        if i > burnin:
            yield y, u
        
    # This can't be right because it doesn't have the numbers 3 or 4 anywhere in it
    # and indeed the plots come out totally fucked
    
    
    
def direct(n=100000, thin=50):
    """
    sample 
      u ~ N(3, 16)
    y|u ~ N(u, 25)
    
    directly from the definitions
    """
    for i in range(n):
        u = 3 + 4*normal()
        y = u + 5*normal()
        yield y, u
        
    # to test u and y
    # (note: S[,1] = u, S[,0] = y)
    
    # you can test u just by looking at it 
    # you can approximately test
    #  y | u
    # by looking at y in small slices by u
    # take var and mean and hists of these to see
    # --> mean should be around the mean of the lower and upper bounds
    # --> var should be 25
    # like so:
    # S[(2 < S[,1]) & (S[,1] < 2.3), 0])
    #
    # and ditto to test u | y
    
    # results:
    # u is correct              (mean = 3, var = 16)
    # y is correct              (mean = 3, var = 41)
    # y | u seems to be correct (mean = u, var = 25)
    # u | y seems to be wrong   (mean = is not directly u, but it is a a linear function of u, with slope around 1/3rd, var = 9.5)
    # so my derivation is fucking wrong