"""
univariate prototyping of flipping conditioning
 in order to be able to do the multivariate gibbs sampler

"""

from numpy import array
from numpy.random import normal

def gibbs1(burnin=10000, n=1000000, thin=50):
    """
    sample 
      u ~ N(3, 16)
    y|u ~ N(u, 25)
    this is a /random effects/ model
    
    via gibbs sampling
    """
    
    u, y = 0,0
    
    for i in range(burnin+n):
        for _ in range(thin):
            y = u + 5*normal()
            u = y - 5*normal() #<-- there's a bug here. the bug is assuming that y and the second normal should be independent.
        if i >= burnin:
            yield u, y
        
    # This can't be right because it doesn't have the numbers 3 or 16 anywhere in it
    # and indeed the plots come out totally fucked
    # 
    
    
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
        yield u, y
        
    # to test u and y
    # (note: S[,0] = u, S[,0] = 1)
    
    # you can test u just by looking at it 
    # you can approximately test
    #  y | u
    # by looking at y in small slices by u
    # take var and mean and hists of these to see
    # --> mean should be around the mean of the lower and upper bounds
    # --> var should be 25
    # like so:
    # S[(2 < S[,0]) & (S[,0] < 2.3), 1])
    #
    # and ditto to test u | y
    
    

def gibbs2(burnin=10000, n=1000000, thin=50):
    """
    sample 
      u ~ N(3, 16)
    y|u ~ N(u, 25)
    
    u | y ~ N
    
    this is a /random effects/ model
    
    via gibbs sampling
    """
    
    u, y = 0,0
    
    for i in range(burnin+n):
        for _ in range(thin):
            y = u + 5*normal()
            u = (5**2*3 + 4**2*y)/(5**2+4**2) + ((4**2 * 5**2)/(4**2 + 5**2))**.5 * normal()
        if i >= burnin:
            yield u, y
        
    # This can't be right because it doesn't have the numbers 3 or 16 anywhere in it
    # and indeed the plots come out totally fucked
    # 
    

def analyze(S):
    """
    given a matrix of joint (u, y) samples (shape: n x 2),
    print means and variances for u, y, u | y and y | u
    the conditionals give a whole table of results--only approximate because you'll never have
    an infinitely thin binning that has any samples in it at all--so you can compare the independent
    variable (i.e. the condition) on the left against the mean and variance for the dependent variable.
    """
    S = array(list(S))
    
    u = S[:,0]
    u = u.mean(), u.var()
    
    y = S[:,1]
    y = y.mean(), y.var()

    y_given_u = [(B, S[(B < S[:,0]) & (S[:,0]<B+1), 1]) for B in range(-100, 100)]
    y_given_u = [(B, a) for B, a in y_given_u if len(a)]
    y_given_u = [(B, a.mean(), a.var()) for B, a in y_given_u]
    y_given_u = array(y_given_u)
    
    u_given_y = [(B, S[(B < S[:,1]) & (S[:,1]<(B+1)), 0]) for B in range(-100, 100)]
    u_given_y = [(B, a) for B, a in u_given_y if len(a)]
    u_given_y = [(B, a.mean(), a.var()) for B, a in u_given_y]
    u_given_y = array(u_given_y)
    
    print("Marginals:")
    
    print("u:")
    print(u)
    print()
    
    print("y:")
    print(y)
    print()
    
    
    print("(approximate) Conditionals:")
    print("y | u:")
    print(y_given_u)
    print()
    
    print("u | y:")
    print(u_given_y)
    print()
    
    return S

    # results:
    # u is correct              (mean = 3, var = 16)
    # y is correct              (mean = 3, var = 41)
    # y | u seems to be correct (mean = u, var = 25)
    # u | y seems to be wrong   (mean ~= u/3 + 2, var = 9.5)
# so my derivation is fucking wrong
# The bug:#
#
# u = 3 + 4*z1
# y = u + 5*z2
# z1 indep z2
#
## I worked out u = y - 5*z2
##  but y and z2 are *not* independent
#  i think if i condition u | y then i also need to condition z2 | y
# . .... but this is.. hard??
# z2 = (y - u)/5
# 
# soooo
#  uh
# u = y - 5(y - u)/5
#   = y - y + u
#   = u
# durp
# that doesn't help
# 


print("gibbs1:")
Sg1 = analyze(gibbs1())

print("Direct:")
Sd = analyze(direct())

# Here's the bug:

# So let's do it from first principles: from the p.d.f.:
# First, definitions, since i didn't actually write these into computor:
# z1, z2 ~iid~ N(0, 1)
# u = 3 + 4*z1
# y = u + 5*z2
# these are chosen so that
#      u ~ N(3, 4**2)
#    y|u ~ N(u, 5**2)
# which is hopefully directly obvious from the definitions.
# So, you'd think then that you could flip it like this:
# u = y - 5*z2
# and therefore to sample u | y, take another N(y, 5**2)
#  but that's absurd because it means that the gibbs sampler totally ignores the 3 and 4 from the marginalized u dist,
#   and analyze(gibbs1()) demonstrates how the doing this makes all the distributions absurd.
# So what's wrong?
# Here's it:
#  y and z2 are not independent, so conditioning on y means screws up z2 as well. I *think* (but do not have a proof handy) that if we knew z2 | y we could, however trying to derive that is an infinite-regress rabbit hole (TODO: explore this experimentally; write a direct2() that records z1 and z2 in addition to y and u and see if there's anything sensible in the conditionals)
# aside: why does y | u work? --> because u only depends on z1, so u *is* independent of z2 and conditioning y | u = u|u + 5*z2|u makes sense.
# 
# So, now that the bug is identified, how do we fix it?
# I wasted a lot of paper showing that u = (u + 5z) - 5z = u over and over again.
# Eventually I decided to try going lower: going to the pdfs.
#
# So let's do that now, being more careful to distinguish the random variables from the values we process them at: capitals are random vars, small letters are constants.
# p(u, y) = 
#   p(U = u, Y = y)
# = p(3 + 4*Z1 = u, u + 5*Z2 = y)  
# = p(Z1 = (u-3)/4, Z2 = (y-u)/5)   #write in standard normal form, which I can look up the pdf for
# = p(Z1 = (u-3)/4) * P(Z2 = (y-u)/5)  #since z1 and z2 are independent
# = C * exp(-1/2 * [(u-3)/4]**2) * exp(-1/2 * [(y-u)/5]**2)  #where C is the normalizing constant ((in this case C is the product of the two normalizing constants for Z1 and Z2 together)
# = C * exp(-{[(u-3)/4]**2 + [(y-u)/5]**2}/2)
# 
# To do gibbs sampling, we are interested in
#   P(y | u)  and  P(u | y)
# which are actually
#   P(U=u, Y=y)/P(U=u)              ( and  )
# = P(U=u, Y=y)/integrate(P(U=u, Y=y)dy)   #by marginalization
# = C * exp(-{[(u-3)/4]**2 + [(y-u)/5]**2}/2)  /  integrate( C * exp(-{[(u-3)/4]**2 + [(y-u)/5]**2}/2) dy )
# = exp(-{[(u-3)/4]**2 + [(y-u)/5]**2}/2)  /  (  exp(-{[(u-3)/4]**2}/2) * integrate( exp(-{[(y-u)/5]**2}/2) dy )
# = exp(-{[(y-u)/5]**2}/2)  /  ( integrate( exp(-{[(y-u)/5]**2}/2) dy )
# here is where my human pattern matching comes in: I recognize the numerator as part of the pdf of a N(u, 5**2),
#   and the denominator as the same but summed--i.e. it is the normalizing constant to match the numerator.
# So this is the pdf of a N(u, 5**2)
# so Y | U=u ~ N(u, 5**2)
# 
# U | Y=y is trickier because the term we need to integrate out--u--appears twice.
# Making it appear only once will make things much easier:
# subproblem: simplify (-with-respect-to-u)
#   [(u-3)/4]**2 + [(y-u)/5]**2
# = (u**2 - 2*3*u + 3**2)/4**2  +  (y**2 - 2*y*u + u**2)/5**2
# = [(5**2 u**2) - (5**2 2*3*u) + (5**2 3**2) + (4**2 y**2) - (4**2 2*y*u) + (4**2 u**2)] / [4**2 5**2]
# = [(4**2 + 5**2)u**2 - 2*(5**2 * 3 + 4**2 * y)u + (5**2 3**2 + 4**2 y**2)]/[4**2 5**2]  #collect, and in a form that is close to a completed square
# = {(4**2 + 5**2)[u**2 - 2*(5**2 * 3 + 4**2 * y)/(4**2 + 5**2)u] + (5**2 3**2 + 4**2 y**2)]}/{4**2 5**2}
# now, complete the sqaure. Since this is large, I'll rederive completing the square and then sub:
#   u**2 - 2*b*u
# = u**2 - 2*b*u + b**2 - b**2 #make clever substitution
# = (u-b)**2 - b**2
# And in this case, b = (5**2 * 3 + 4**2 * y)/(4**2 + 5**2), a term which *only depends on y* (so that the integral can treat it as a constant)
# = {(4**2 + 5**2)[(u-b)**2 - b**2] + (5**2 3**2 + 4**2 y**2)]}/{4**2 5**2}
# = {(4**2 + 5**2)(u-b)**2 - (4**2 + 5**2)b**2 + (5**2 3**2 + 4**2 y**2)]}/{4**2 5**2}
# = (4**2 + 5**2)(u-b)**2/{4**2 5**2} - (4**2 + 5**2)b**2/{4**2 5**2} + (5**2 3**2 + 4**2 y**2)]/{4**2 5**2}   #split the terms because holla
# = (u-b)**2/({4**2 5**2}/(4**2 + 5**2)) - (4**2 + 5**2)b**2/{4**2 5**2} + (5**2 3**2 + 4**2 y**2)]/{4**2 5**2}   # tweak the factors to collect the variance term into one unit
# = 
#  ^ this barely looked any better in handwriting.

# actually, I can handwave the rest of the derivation. There's a zillion little fiddly terms and none of them will end up mattering because pdf has to be normalized to sum to 1, so all those non-u terms will just end up as part of the normalization.
# Deriving P(u | y) from here will be identical to doing P(y | u), but with which terms are considered constant and which are variable tossed up.
#  so we'll end up with some normal distribution again, and we just need to figure out which one.
# In the above formula we have all the information that uniquely determines which normal dist it is, because the e^(-1/2) part is fixed for all normals. It's what's inside of that function that is important.
# the mean is the amount that u is shifted
# the variance is the amount the shifted u is scaled (in full, not in square root--that would be the std.dev.)
# namely:
# mean(U | Y=y) = b = (5**2 * 3 + 4**2 * y)/(4**2 + 5**2)
#  var(U | Y=y) =               {4**2 5**2}/(4**2 + 5**2)

# so we should have
# y | u ~ N(u, 5**2), as before, but ALSO, that
# u | y ~ N((5**2*3 + 4**2*y)/(5**2+4**2),  (4**2 * 5**2)/(4**2 + 5**2))
# AHA
# AWESOME
# THIS MATCHES THE DATA. Try it out:
# the variance for u | y is about 9.75, very close to most of the empirical conditional variances,
# and the mean, as a function of y, has a slope of about .4 and intercept of about 2

print("gibbs2:")
Sg2 = analyze(gibbs2())

