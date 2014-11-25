Questions
=========

* [ ] What parts of the assignment are typos?
* [ ] Why is the U | Y formula in the assignment backwards from our derivation?
  * Weirdly, it seems to check out: `gibbs_hierarchical_normal_unrolled()` (ours) gives empirically the same distribution as `gibbs_hierarchical_normal_wrong()` (lysy's).
* [ ] How do you actually *use* this? Sure, I can sample a bunch of conditionals, but how does that help in estimation(/regression/inference/whatever)? Is there a secret bayesian posterior to sample that I'm missing, here?
  * I doubt the method is: pick a B, sample, compare results, output B that gives the least error, because there's infinite Bs to check.
  * Probably reading the Everson & Morris (2002) paper will shed light on this.
* [ ] Why bother to gibbs-sample at all if we know the *full* covariance matrix? We have a much faster and parallelizable method for multinormals!
  * This is especially confusing since I don't see any way to work out the conditional formulas without incidentally working out the full covariance matrix.
