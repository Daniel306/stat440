440 Ideas:

MultiD Distribution Comparison
==============================

In 1d, there's qqplots and the k-s test,
but how do we do something like eg test
that the InvWishart code is behaving?

 -> plot a scatterplot matrix?
 -> plot ????????????????


Corpus Management by Machine Learning
=====================================
 - a tool for dealing with large heterogenous corpii
  namely, pirated media
   - music (taking music similarities)
   - epubs (taking text similarities)
   - 


InvWishart Project
==================
(via mlysy)

Implement the Inverse Wishart distribution for R and Scipy. Use Rcpp for R and Pyrex for Python.



Tumblr
=======

Classifier to discriminate posts a person wrote from reblogs. It is a bit tricky because reblogs can include comments, and those can be mini-essays in themselves.
But sometimes (often) the reblogs are shit.
But sometimes maybe you want to look at only the shit.
THE POSSIBILITIES ARE ENDLESS.

Use case:
 > extract everywhere someone is gossiping about you
 > filter out fandom crud (without relying on the user to tag their posts correctly)

Confidence Interval Project
===========================

How do we compare confidence intervals? A CI (or the bayesian analogue) gives an region (and interval) which has high probability of covering the quantity in question.
So, 

See: Terry's brain model talk.


Multinomial Project
===================
 - defining categorical data
   - if you have tagged data, that is, each of your n data points comes with a set {t1, t2, t3, ..} \subset {t1, .. t_j .. tK}
     there are two ways to approach it:
       - treat the tags as K distinct multinomial (actually, binomial) problems (so, there are K {iid Ber(p_j)} dist'ns, or K ind Bin(n, p_j) dist'ns)
       - treat the Bernoullis as a binary string (that is, {t1,t4,t6, ...} = {yes, no, no, yes, no, yes, ....}) and treat those binary strings as multinomial (from a 2^k-sized multinomial distribution)

      the latter is more tractable: it's easy to do the mapping, easy to estimate frequencies of the multi, but loses a lot of information in the process; 
        because of the curse of dimensionality it is typical that a sample with n >> k (so n/k >> 1) 
        becomes n/k << 1    
   
   How do we define correlation between tags? --> tags are uncorrelated when P(t1, t2) = P(t1)P(t2). So, you can use bayes rule:
     P(t1|t2) = P(t1, t2)/P(t2).
     If P(t1|t2) = P(t1) then t1 is uncorrelated with
     (and this implies t2 is uncorrelated with t1 because
         P(t2|t1) = P(t1|t2)P(t2)/P(t1) = P(t1)P(t2)/P(t1) = P(t2))
     ( interestingly, we don't have reflexivity because P(t1|t1) = P(t1, t1)/P(t1) = 1
          intuitively, correlated means "as you change t2 you change t1" and uncorrelated is "as you change t2 you don't change t1" 
           and in P(t1|t1) = 1 we do have that: as you change t1 it's still 1
         ---huh, does this contradict the other def'n of uncorrelated (which is cov(x,y)=0) or the def'n of indepedent (which is p(a,b) = p(a)p(b) 
 
     P(t1,t1) != P(t1)P(t1) --- the right side is smaller (unless P(t1)=1)
     )

   --> some work to relax the condition from strict independence to near independence
     and to 

   Is there something in between?
        
 - get corpii by sampling:
   - pixiv
   - wikipedia
   - our own music libraries
   - soundcloud
   - tumblr (hashtags)
   - twitter (hashtags)
   - facebook (hashtags)
   - 
 -
