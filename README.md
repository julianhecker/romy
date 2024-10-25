
<!-- README.md is generated from README.Rmd. Please edit that file -->



# ROMY: Robust Omics MethodologY

<!-- badges: start -->
<!-- badges: end -->

The *romy* R package provides flexible and robust tools for the analysis of multi-omics data. **This package is still under development and not for use yet.**

## Installation

You can install *romy* using:

``` {.r}
#library(devtools)
#install_github("julianhecker/romy")
library(romy)
```

## Background

The *romy* R package implements five different analysis components for the analysis of multi-omics data.
- association testing: testing for association between different omics factors and phenotype variables while adjusting for covariates/confounders Z.
- quantile association testing: testing for association between different omics factors and phenotype variables while adjusting for covariates/confounders Z, using an quantile-based approach. Inspired by optimal transport. - (co)variance testing: testing for an effect of X on the covariance between two factors in Y, or the variance of one factor in Y, adjusting for Z.
- component testing: testing groups of variables in X for association with Y, adjusting for Z.
- interaction testing: testing for interaction between factor in X and factor in Z in their effect on Y.

All these approaches utilize the concepts of debiased machine learning and/or robust statistics (Chernozhukov et al. 2018, Shah and Peters 2020).

## Component 1: association testing

The first analysis approach in *romy* performs association testing between variables in Y and X, adjusting for Z. The measurements in X and Y can be transformed using data transformation functions. In this case, all pairs of transformations are tested separately and jointly. The separate p-values are also combined in an ACAT (Aggregated Cauchy Association Test).

``` {.r}
X=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Z=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Y=matrix(rnorm(1000*30), nrow=1000, ncol=30)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

at=association_testing(Y, X, Z)
head(at)
#>   X_id Y_id      p_value
#> 1   X1   Y1 6.337642e-01
#> 2   X2   Y1 4.677403e-01
#> 3   X3   Y1 8.090326e-01
#> 4   X1   Y2 8.770232e-05
#> 5   X2   Y2 5.985061e-01
#> 6   X3   Y2 9.537739e-02
```

## Component 2: association testing based on quantiles

The second analysis approach in *romy* performs association testing between variables in Y and X, adjusting for Z, and using quantile regression. We implemented a test statistic inspired by optimal transport. The approach performs three different tests using three different quantile levels. There is also a joint test and an ACAT based on the three separate p-values.

``` {.r}
X=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Z=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Y=matrix(rnorm(1000*30), nrow=1000, ncol=30)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

at=association_testing_quantile(Y, X, Z)
head(at)
#>   X_id Y_id p_values1 p_values2 p_values3 p_overall    p_ACAT
#> 1   X1   Y1 0.6084152 0.8510714 0.9199801 0.9386067 0.8921246
#> 2   X2   Y1 0.3278302 0.2590433 0.7780105 0.5118511 0.4750165
#> 3   X3   Y1 0.9886504 0.6938433 0.3791579 0.8509580 0.9582336
#> 4   X1   Y2 0.4366106 0.5139260 0.9204946 0.6291262 0.7574386
#> 5   X2   Y2 0.9715949 0.3598700 0.1980188 0.5478770 0.8726927
#> 6   X3   Y2 0.3766123 0.9378681 0.3845694 0.6441709 0.7774437
```

## Component 3: (co)variance testing

The third analysis approach aims to test for effects of variables in X on the (co)variance in Y, adjusting for Z.

``` {.r}
X=matrix(rnorm(1000*3),nrow=1000,ncol=3)
Z=matrix(rnorm(1000*3),nrow=1000,ncol=3)
Y=matrix(rnorm(1000*10),nrow=1000,ncol=10)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

ct=covar_testing(Y, X, Z)
head(ct)
#>   X_id Y_id1 Y_id2    p_value
#> 1   X1    Y1    Y1 0.01260383
#> 2   X2    Y1    Y1 0.25025525
#> 3   X3    Y1    Y1 0.83762208
#> 4   X1    Y1    Y2 0.64211310
#> 5   X2    Y1    Y2 0.63997879
#> 6   X3    Y1    Y2 0.14148586
```

## Component 4: component testing

The fourth analysis approach aims to test structured groups in X (sets of variables) for association with Y, adjusting for Z. This utilizes the first analysis approach and a factor analysis to reduce dimensionality.

``` {.r}
X=matrix(rnorm(1000*30),nrow=1000,ncol=30)
Z=matrix(rnorm(1000*5),nrow=1000,ncol=5)
Y=matrix(rnorm(1000),nrow=1000,ncol=1)
colnames(Y)=paste0("Y",1:ncol(Y))
colnames(Z)=paste0("Z",1:ncol(Z))

#ct=component_testing(Y, X, Z)
```

## Component 5: interaction testing

The fifth analysis approach aims to perform robust interaction testing between variables in X and Z in their effect on Y. This is similar to the RITSS (Robust Interaction Testing using Sample Splitting) method.

``` {.r}
X=matrix(rnorm(1000*30),nrow=1000,ncol=30)
Z=matrix(rnorm(1000*5),nrow=1000,ncol=5)
Y=matrix(rnorm(1000*1),nrow=1000,ncol=1)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))
colnames(Z)=paste0("Z",1:ncol(Z))

it=interaction_testing(Y, X, Z)
head(it)
#>   X_id Y_id    p_value
#> 1   X1   Z1 0.98686588
#> 2   X2   Z1 0.05824525
#> 3   X3   Z1 0.12079220
#> 4   X4   Z1 0.66926188
#> 5   X5   Z1 0.19661483
#> 6   X6   Z1 0.43007784
```

## References

Bai and Ng 2002. Determining the Number of Factors in Approximate Factor Models.
Li and Li 2002. Integrative Factor Regression and Its Inference for Multimodal Data Analysis.
Hecker et al. 2022. A robust and adaptive framework for interaction testing in quantitative traits between multiple genetic loci and exposure variables.
Shah and Peters 2020. The Hardness of Conditional Independence Testing and the Generalised Covariance Measure.
Chernozhukov et al. 2018. Double/debiased machine learning for treatment and structural parameters.
