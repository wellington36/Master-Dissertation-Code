## Code for the SBC in Conway-Maxwell-Poisson Models chapter

### Replicate Figures
Dependencies: `com_poisson_pmf.R` and `rcomp_benson.R`.
- Fig. 6: `mean_and_variance_check.R`
- Fig. 7: `marginal_joint_proposed_prior.R`
- Table 2: `packagers_vs_our_method.R`

### Other
- `check_hyperparameters.R`: plot distribution of CMP given the hyperparameters for gamma prior distributions.
- `rcomp.R`: generate data from CMP based in Flavio's proposed idea.
- `model_check_monte_carlo.R`: run Monte Carlo method to evaluate the mean and variance with the total expectation and total variance law.
