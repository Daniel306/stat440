Questions
=========

[ ] Explain the typos in the assignment.
[ ] Explain why the conditional formula as given works. It is backwards from what I derived, but empirically it seems to check out.
[ ] How do you actually *use* this? Sure, I can sample a bunch of conditionals, but how does that help in estimation(/regression/inference/whatever)? Is there a secret bayesian posterior to sample that I'm missing, here? I doubt the method is: pick a B, sample, compare results, output B that gives the least error, because there's infinite Bs to check.
[ ] Why bother to gibbs-sample if the same formulas we use to derive the component matrices also give us the *full* covariance matrix? We have a much faster and parallelizable method! Just sample the multivariate normal on the larger covariance matrix, then partition it.