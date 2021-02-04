# Sequential discovery models (class `sdm`)

Sequential discovery models are a Bayesian method to construct, fit and predict the accumulation curves arising from the frequencies with which species are observed. For a theoretical description, see:

 * Zito, A., Rigon, T., Ovaskainen, O. and Dunson, D. B. (2020+): [Bayesian nonparametric modelling of sequential discoveries](https://arxiv.org/abs/2011.06629)
 
The models available for the latent variables are the three-parameter log-logistic distribution (`"LL3"`, *the default*) and the Weibull distribution (`"Weibull"`). Both methods assume that the *asymptotic species richness*, which is the total number of species observable in the sample, is always finite. Such an assumption allows to determine how close the accumulation curves are to convergence according to the selected model. **Note**: the models provide reliable results when the frequencies are 

As a working example, we consider the 

