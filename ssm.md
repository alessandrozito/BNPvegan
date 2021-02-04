# Species sampling models (`ssm` class)

Species sampling models are a Bayesian nonparametric tools that have a large variety of applications in ecology. For a general overview, one can refer to:

* Lijoi, A., Mena, R.H., Pruenster, I. (2007) [Bayesian nonparametric estimation of the probability of discovering new species](https://academic.oup.com/biomet/article-abstract/94/4/769/246082), *Biometrika* **94**(4), 769--786.

* De Blasi, P., Favaro, S., Lijoi, A., Mena, R.H., Pruenster, I., Ruggiero, M.: (2015): [Are Gibbs-type priors the most natural generalization of the Dirichlet process?](https://arxiv.org/abs/1503.00163), *IEEE Transactions on Pattern Analysis and Machine Intelligence* **37**(2), 212-229.

As a working example, we make use of the frequencies available in the `data("Lepidoptera")` dataset. The sampling scheme currently availabe are the Dirichlet process `DP` and the Pitman-Yor process `PY`.

```r 
# Load the both the vegan and BNP vegan libraries into memory
library(vegan) 
library(BNPvegan)
data("Lepidoptera")
```

The `Lepidoptera` dataset contains a vector of frequencies, denoting the number of times a given species has been observed.

```r 
frequencies <- as.numeric(Lepidoptera) # Frequencies of each species
```
The whole dataset is displayed in the following for the sake of completeness

```r 
 [1]   1   4  28 604   5  32   4   8   2  13   4  32 190   1  51  53   1
 [18] 221  22   6   3  57   3   6   1 246   4  18  53  43  23  49  19  34
 [35]  10  22  23  20   1   1  64  16   8  60   2  32  16 743  16  73   6
 [52]  11  10   8 604  34   7   3  23   6   1 129 181   5   1   7   3   5
 [69]   3  61 743  14   8   6   1 221   6   4   2  67  12  25   4   1  20
 [86]  19   1  99   6   3  15   6  22  44 190   2   1  15   1 743  52   4
[103]   8  28   5   1   1  44   1  44 187  44 109   1   4   3  36  58  21
[120]   6   4   3  49   3  51  12   4 154  44  49   1  24  44  28   1   5
[137]   5   7 109   4  15   1   1  13  51 112   5  23  49  16   4   2  22
[154]   7   4   1   2   4   6 187   1  57   1  13  19   2   7  34  21   5
[171] 306   4   5  22  22   1   2   7   1   6  44 181   4  44   7  34  53
[188] 148   8 109   9  28  19   2   3   3   4   5   6 572   3  24  32   1
[205]  28 154  54  53  11  99   3  20   4  64  28   1  10  12   4   4  96
[222] 221  96 333   1   5 112   1   1   9  25  27 464   6  18  43 112  10
[239] 154 246
```

## Estimating the sample coverage

## Bio-diversity

## Estimating the rarefaction

## Extrapolation of the curve
