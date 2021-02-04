# Sequential discovery models (class `sdm`) - :warning: Under development

Sequential discovery models are a Bayesian method to construct, fit and predict the accumulation curves arising from the frequencies with which species are observed. For a theoretical description, see:

 * Zito, A., Rigon, T., Ovaskainen, O. and Dunson, D. B. (2020+): [Bayesian nonparametric modelling of sequential discoveries](https://arxiv.org/abs/2011.06629)
 
The models available for the latent variables are the three-parameter log-logistic distribution (`"LL3"`, *the default*) and the Weibull distribution (`"Weibull"`). Both methods assume that the *asymptotic species richness*, which is the total number of species observable in the sample, is always finite. Such an assumption allows to determine how close the accumulation curves are to convergence according to the selected model. **Note**: the models provide reliable results when the frequencies are 

As a working example, we consider the following examples of accumulation curve generated from a set of frequencies of fungal operational taxonomic units (OTU).

```R
# Load the frequencies
frequencies <- fungalOTU

# Fit the model with 1000 resamples
set.seed(1)
fit <- sdm(frequencies, n_resamples = 1000, model = "LL3")
```

To summarize the output, just run
```R
summary(fit)
```
The result is
```R
Model:
	 Three-parameter log-logistic (LL3)
	 Number of resamples: 1000

Quantities:
	 Abundance: 21243
	 Richness: 563
	 Expected species at infinity: 683
	 Standard deviation at infinity: 26.12
	 Sample saturation: 0.8243

Parameters:
	     alpha       sigma         phi     loglik
	 ---------  ----------  ----------  ---------
	  91.86718   0.0269072   0.9999871   -1975.74
```
To make a plot of the chosen accumulation curve, run
```R
plot(fit)
```


Finally, to extract relevant quantities, do:
```R
coef(fit)
logLik(fit)
asymptotic_richness(fit)
```

