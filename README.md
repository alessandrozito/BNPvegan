# The `BNPvegan` R package - :warning: under development

This **R** package is an implementation of various Bayesian nonparametric models for the analysis of ecological data. The package consists of two main methods: **species sampling models** (`ssm`) and  **species discovery models** (`sdm`).



The `BNPvegan` package can be installed by running the following commands:

```r 
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("alessandrozito/BNPvegan")
```


## Species sampling models (`ssm` class)

Species sampling models are a Bayesian nonparametric tools that have a large variety of applications in ecology. For a thorough description of the main models, see:

* De Blasi, P., Favaro, S., Lijoi, A., Mena, R.H., Pruenster, I., Ruggiero, M.: (2015): [Are Gibbs-type priors the most natural generalization of the Dirichlet process?](https://arxiv.org/abs/1503.00163), *IEEE Transactions on Pattern Analysis and Machine Intelligence* **37**(2), 212-229.


As a working example, we make use of the frequencies available in the `data("Lepidoptera")` example. The sampling scheme availabe are the Dirichlet process `DP` and the Pitman-Yor process `PY`

```r 
# Load the library into memory
library(BNPvegan)
data("Lepidoptera")
```

## Sequential discoveries models (`sdm` class) 

Sequential discovery models are a Bayesian method to construct, fit and predict the accumulation curves arising from the frequencies with which species are observed. For a theoretical description, see:

 * Zito, A., Rigon, T., Ovaskainen, O. and Dunson, D. B. (2020+): [Bayesian nonparametric modelling of sequential discoveries](https://arxiv.org/abs/2011.06629)
 
 The models ava

```r 
# Add here relevant commands
```
