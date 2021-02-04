# The `BNPvegan` R package - :warning: Under development

This **R** package is an implementation of various Bayesian nonparametric models for the analysis of ecological data. The **R** package is currently under development. The package consists of two main methods: **species sampling models** (`ssm`) and  **species discovery models** (`sdm`).

The `BNPvegan` package can be installed by running the following commands:

```r 
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("alessandrozito/BNPvegan")
```

We discuss the two class of models and the associated functions in the following tutorials. Please note that there might be some overlap between these two modeling classes. 

- In the [**SSM tutorial**](ssm.md) we discuss **species sampling models** and their usage in ecology. Further documentation is also given in these [**slides**](slides_SSM.pdf).

- In the [**SDM tutorial**](sdm.md) tutorial we discuss **species discovery models** and their usage in ecology.
