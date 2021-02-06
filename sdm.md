# Sequential discovery models (class `sdm`) - :warning: Under development

Sequential discovery models are a Bayesian method to construct, fit and predict the accumulation curves arising from the frequencies with which species are observed. For a theoretical description, see:

 * Zito, A., Rigon, T., Ovaskainen, O. and Dunson, D. B. (2020+): [Bayesian nonparametric modelling of sequential discoveries](https://arxiv.org/abs/2011.06629)
 
The models available for the latent variables are the three-parameter log-logistic distribution (`"LL3"`, *the default*) and the Weibull distribution (`"Weibull"`). Both methods assume that the *asymptotic species richness*, which is the total number of species observable in the sample, is always finite. Such an assumption allows to determine how close the accumulation curves are to convergence according to the selected model. **Note**: the models provide reliable results when the sample size is large enough (i.e. larger  than 5000). 

As a working example, we consider the following examples of accumulation curve generated from a set of frequencies of fungal operational taxonomic units (OTU), called `fungalOTU`. Notice that the construction of the accumulation curves is inherently order dependent. To cope for this fact, we adopt a resampling approach by specifying the parameter `n_resamples = 1000`.

 This ensures a reasonable computational time irrespective of the size of the curve. At every resample, the model samples one random sequence of discoveries from the observed `frequencies` and runs the `model` specified. The curve and the parameters returned for which the chosen `model` are the ones corresponding to the median asymptotic species richness across reamples. To ensure exact reproducibility, it is recommended to set a seed before running (under sufficiently `n`, differences in saturation are minimal). 

```R
library(BNPvegan)

# Load the frequencies
frequencies <- fungalOTU

# Fit the model with 1000 resamples
set.seed(1) 
fit <- sdm(frequencies, model = "LL3", verbose = TRUE, n_resamples = 1000)
```

To summarize the output, just run
```R
summary(fit)

Model:
	 Three-parameter log-logistic (LL3)
	 Number of resamples: 500

Quantities:
	 Abundance: 21243
	 Richness: 563
	 Expected species at infinity: 689
	 Standard deviation at infinity: 26.24
	 Expected new species to discover: 126
	 Sample saturation: 0.8171

Parameters:
	     alpha       sigma         phi     loglik
	 ---------  ----------  ----------  ---------
	  40.52234   0.1385731   0.9999812   -2081.25
```

To make a plot of the chosen accumulation curve, run:
```R
plot(fit, type = "rarefaction")  # Note: default type is rarefaction, which plots also the observed accumulation curve
```

<img src="https://github.com/alessandrozito/BNPvegan/blob/master/img/sdm_plot.png" width="600" >

To plot the out of sample prediction, just type
```R
n <- length(fit$discoveries)
plot(fit, type = "extrapolation", m = n)  # m are the additionall points to compute the prediciton. Default is m=n
```

<img src="https://github.com/alessandrozito/BNPvegan/blob/master/img/sdm_plot_extrapolation.png" width="600" >

To extract relevant quantities, do:
```R
coef(fit)                   # Coefficients alpha, sigma and phi
logLik(fit)                 # Loglikelihood of the best curve
asymptotic_richness(fit)    # asymptotic species richness and saturation
```

Finally, rarefactions and predictions can be obtained as follows 
```R
# To obtain a generic prediction with the estimated parameters 
newdata <- c(10, 100, 1000)
predict(fit, newdata = newdata)

# To obtain the rarefaction curve, do one of the following (they are equivalent. predict is a generic prediction method)
predict(fit)
rarefaction(fit)

# To obtain the extrapolation, use
m = 10
extrapolation(fit, m = m)  # extrapolates the curve for additional 10 point.
```





