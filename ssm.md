## Species sampling models (`ssm` class)

Species sampling models are a Bayesian nonparametric tools that have a large variety of applications in ecology. For a general overview, one can refer to:

* Lijoi, A., Mena, R.H., Pruenster, I. (2007) [Bayesian nonparametric estimation of the probability of discovering new species](https://academic.oup.com/biomet/article-abstract/94/4/769/246082), *Biometrika* **94**(4), 769--786.

* De Blasi, P., Favaro, S., Lijoi, A., Mena, R.H., Pruenster, I., Ruggiero, M.: (2015): [Are Gibbs-type priors the most natural generalization of the Dirichlet process?](https://arxiv.org/abs/1503.00163), *IEEE Transactions on Pattern Analysis and Machine Intelligence* **37**(2), 212-229.

As a working example, we make use of the frequencies available in the `data("Lepidoptera")` dataset. The sampling scheme currently availabe are the Dirichlet process `DP` and the Pitman-Yor process `PY`.

```r 
# Load the library into memory
library(BNPvegan)
data("Lepidoptera")
```

The `Lepidoptera` dataset contains a vector of frequencies, denoting the number of times a given species has been observed. Hence


