# Sequential discovery models (class `sdm`) - :warning: Under development

Sequential discovery models are a Bayesian method to construct, fit and predict the accumulation curves arising from the frequencies with which species are observed. For a theoretical description, see:

 * Zito, A., Rigon, T., Ovaskainen, O. and Dunson, D. B. (2020+): [Bayesian nonparametric modelling of sequential discoveries](https://arxiv.org/abs/2011.06629)
 
The models available for the latent variables are the three-parameter log-logistic distribution (`"LL3"`, *the default*) and the Weibull distribution (`"Weibull"`). Both methods assume that the *asymptotic species richness*, which is the total number of species observable in the sample, is always finite. Such an assumption allows to determine how close the accumulation curves are to convergence according to the selected model. **Note**: the models provide reliable results when the sample size is larger than 5000. 

As a working example, we consider the following examples of accumulation curve generated from a set of frequencies of fungal operational taxonomic units (OTU), called `fungalOTU`. Notice that the construction of the accumulation curves is inherently order dependent. To cope for this fact, we adopt a resampling approach by specifying the parameter `n_resamples`. The default value is set to `n_resamples=1000`, but should be reduced under large datasets to shorten the computational time. At every iteration, the model samples one random sequence of discoveries from the observed `frequencies` and runs the `model` specified. The output corresponds to the curve and the parameters for which the chosen `model` reaches the highest log-likelihood. To ensure exact reproducinbility, it is recommended to set a seed before running. 

```R
# Load the frequencies
frequencies <- fungalOTU

# Fit the model with 1000 resamples
set.seed(1) 
fit <- sdm(frequencies, n_resamples = 1000, model = "LL3", verbose = TRUE)
```

To summarize the output, just run
```R
summary(fit)

Model:
	 Three-parameter log-logistic (LL3)
	 Number of resamples: 1000

Quantities:
	 Abundance: 21243
	 Richness: 563
	 Expected species at infinity: 682
	 Standard deviation at infinity: 26.11
	 Sample saturation: 0.8255

Parameters:
	     alpha       sigma         phi     loglik
	 ---------  ----------  ----------  ---------
	  91.47824   0.0275633   0.9999869   -1975.74
```

To make a plot of the chosen accumulation curve, run.
```R
plot(fit)
```

```
![See the plot here](/images/sdm_curve.pdf )
```

Finally, to extract relevant quantities, do:
```R
coef(fit)                   # Coefficients alpha, sigma and phi
logLik(fit)                 # Loglikelihood of the best curve
asymptotic_richness(fit)    # asymptotic species richness and saturation
```

