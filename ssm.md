# Species sampling models (`ssm` class)

Species sampling models are a Bayesian nonparametric tools that have a large variety of applications in ecology. For a **technical** overview, one can refer to:

* Lijoi, A., Mena, R.H., Pruenster, I. (2007) [Bayesian nonparametric estimation of the probability of discovering new species](https://academic.oup.com/biomet/article-abstract/94/4/769/246082), *Biometrika* **94**(4), 769--786.

* De Blasi, P., Favaro, S., Lijoi, A., Mena, R.H., Pruenster, I., Ruggiero, M.: (2015): [Are Gibbs-type priors the most natural generalization of the Dirichlet process?](https://arxiv.org/abs/1503.00163), *IEEE Transactions on Pattern Analysis and Machine Intelligence* **37**(2), 212-229.

As a working example, we make use of the frequencies available in the `data("Lepidoptera")` dataset. The sampling scheme currently availabe are the Dirichlet process `DP` and the Pitman-Yor process `PY`.

```r 
# Load the both the vegan and BNP vegan libraries into memory
library(vegan) 
library(BNPvegan)
data("Lepidoptera")
```

The `Lepidoptera` dataset contains a **vector of frequencies**, denoting the number of times a given species has been observed.

```r 
frequencies <- as.numeric(Lepidoptera) # Frequencies of each species

# Other quantities
n <- sum(frequencies) # Sample size, i.e. the abundance
K <- length(frequencies) # Number of distinct species, i.e. the richness. 
M <- as.numeric(table(factor(frequencies, levels = 1:n))) # Frequency of frequencies
```

In this specific dataset, the sample size equals `n = 12548` whereas the richness is `K = 240`. The vector `M` include the **frequencies of frequencies**. In particular, `M[1]` represents the number of singletons, which in this case equals `31`, namely the number of species that have been observed only once. 

## The `ssm` function

The main function for the estimation of species sampling model is called `ssm`. Although several species sampling models existing, we have currently implemented only the `DP` and the `PY`. Please note that the latter is a generalization of the former. 

In order to estimate, say, a `DP` model you can use the following **R** commands:

```r
fit_DP <- ssm(frequencies, "DP") # An object of class ssm
alpha_hat <- fit_DP$param # Maximum likelihood estimate of the parameter
```

This obtains the maximum likelihood estimate for the parameter od the `DP`. In this case we have that `alpha_hat = 41.99451`. Most of the relevant quantities can be obtained using the `summary` function. 

```r
summary(fit_DP)

Model: Dirichlet Process
	 Abundance: 12548
	 Richness: 240
	 alpha: 41.9945
	 Sample coverage: 0.9967
	 Expected species after additional 12548 samples: 269
	 New expected species after additional 12548 samples: 29
	 Posterior Gini diversity: 0.975
```


## Estimating the sample coverage

## Bio-diversity

## Estimating the rarefaction

## Extrapolation of the curve
