# clonealign 2.0

Clonealign 2.0 contains several updated changes, both to the core model and inference, as well as to 
package functionality.

## Major changes:

* While the modelled expectation remains identical, the likelihood employed is now multinomial
rather than negative binomial. This speeds up inference and removes the need for size factor
and variance calculation.
* The preferred package entrypoint is the run_clonealign function, that will run clonealign across
a range of initial parameter values and return the fit that maximizes the ELBO. 

## New functionality

* Post-hoc calculated correlations between the copy number and gene expression are calculated and
assigned to the $correlations slot

## Changes to the core model

* Multinomial likelihood rather than negative binomial
* mu[1] no longer constrained to be 1 to improve optimization in some cases
* It is no longer necessary to specify size factors as the multinomial distribution
implicitly conditions on the total counts per cell


## Changes to inference

* Convergence is monitored by looking at the average change in the previous 10
iterations rather than the single previous iteration, which can be sensitive to random
fluctuations