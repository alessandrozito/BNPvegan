# Sequential discovery models (class `sdm`) - :warning: Under development

Sequential discovery models are a Bayesian method to construct, fit and predict the accumulation curves arising from the frequencies with which species are observed. For a theoretical description, see:

 * Zito, A., Rigon, T., Ovaskainen, O. and Dunson, D. B. (2020+): [Bayesian nonparametric modelling of sequential discoveries](https://arxiv.org/abs/2011.06629)
 
The models available for the latent variables are the three-parameter log-logistic distribution (`"LL3"`, *the default*) and the Weibull distribution (`"Weibull"`). Both methods assume that the *asymptotic species richness*, which is the total number of species observable in the sample, is always finite. Such an assumption allows to determine how close the accumulation curves are to convergence according to the selected model. **Note**: the models provide reliable results when the sample size is large enough (i.e. larger  than 5000). 

As a working example, we consider the following accumulation curve generated from a set of frequencies of fungal operational taxonomic units (OTU), called `fungalOTU`. As a first step, the function constructs the exact rarefaction obtained from the given `frequencies` according the rarefaction equation on page 7/25 in the explanatory  [**slides**](slides_SSM.pdf). Then, it runs the chosen `model` on such a curve, and returns the asymptotic estimates for the species richness and the sample saturation. 

```R
library(BNPvegan)

# Load the frequencies
frequencies <- fungalOTU

# Fit the model
fit <- sdm(frequencies, model = "LL3", verbose = TRUE)
```

To summarize the output, just run
```R
summary(fit)

Model:
	 Three-parameter log-logistic (LL3)

Quantities:
	 Abundance: 21243
	 Richness: 563
	 Expected species at infinity: 690
	 Standard deviation at infinity: 26.26
	 Expected new species to discover: 127
	 Sample saturation: 0.8159

Parameters:
	     alpha       sigma         phi      logLik
	 ---------  ----------  ----------  ----------
	  45.51308   0.1225816   0.9999824   -2067.835
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





