
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ROMY: Robust Omics MethodologY

<!-- badges: start -->
<!-- badges: end -->

The *romy* R package provides flexible and robust tools for the analysis
of multi-omics data. **!This package is still under development.!**

## Installation

You can install *romy* using:

``` r
#library(devtools)
#install_github("julianhecker/romy")
library(romy)
```

## Background

The *romy* R package implements four different analysis components for
the analysis of multi-omics data.  
- Conditional Independence Testing (CIT): testing for association
between different omics factors and phenotype variables while adjusting
for covariates/confounders Z.  
- (co)variance testing (COV): testing for an effect of X on the
covariance between two factors in Y, or the variance of one factor in Y,
adjusting for Z.  
- interaction testing (INTER): testing for interaction between factor in
X and factor in Z in their effect on Y.  
- Latent Factor (LF): Estimation of latent factors that correspond to
heterogenous effect profiles.

The first three approaches utilize the concepts of debiased machine
learning and/or robust statistics (Chernozhukov et al. 2018, Shah and
Peters 2020).  
The last component (LF) is based on the work of Chen (2022) (under
development).

## CIT

The first analysis approach in *romy* performs association testing
between variables in Y and X, adjusting for Z. The measurements in X and
Y can be transformed using data transformation functions. In this case,
all pairs of transformations are tested separately and jointly. The
separate p-values are also combined in an ACAT (Aggregated Cauchy
Association Test).

``` r
X=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Z=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Y=matrix(rnorm(1000*30), nrow=1000, ncol=30)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

cit=romy_cit(Y, X, Z)
head(cit)
#>   X_id Y_id   p_value
#> 1   X1   Y1 0.2534191
#> 2   X2   Y1 0.9138961
#> 3   X3   Y1 0.2084404
#> 4   X1   Y2 0.4166797
#> 5   X2   Y2 0.2359576
#> 6   X3   Y2 0.2176470
```

## COV

The second analysis approach aims to test for effects of variables in X
on the (co)variance in Y, adjusting for Z.

``` r
X=matrix(rnorm(1000*3),nrow=1000,ncol=3)
Z=matrix(rnorm(1000*3),nrow=1000,ncol=3)
Y=matrix(rnorm(1000*10),nrow=1000,ncol=10)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

cov=romy_cov(Y, X, Z)
head(cov)
#>   X_id Y_id1 Y_id2    p_value
#> 1   X1    Y1    Y1 0.04940465
#> 2   X2    Y1    Y1 0.94595861
#> 3   X3    Y1    Y1 0.88668326
#> 4   X1    Y1    Y2 0.40299703
#> 5   X2    Y1    Y2 0.69135872
#> 6   X3    Y1    Y2 0.46705955
```

## INTER

The third analysis approach aims to perform robust interaction testing
between variables in X and Z in their effect on Y. This is similar to
the RITSS (Robust Interaction Testing using Sample Splitting) method.

``` r
X=matrix(rnorm(1000*30),nrow=1000,ncol=30)
Z=matrix(rnorm(1000*5),nrow=1000,ncol=5)
Y=matrix(rnorm(1000*1),nrow=1000,ncol=1)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))
colnames(Z)=paste0("Z",1:ncol(Z))

inter=romy_inter(Y, X, Z)
head(inter)
#>   X_id Y_id   p_value
#> 1   X1   Z1 0.1135637
#> 2   X2   Z1 0.2610785
#> 3   X3   Z1 0.1089153
#> 4   X4   Z1 0.9121250
#> 5   X5   Z1 0.7975057
#> 6   X6   Z1 0.4250315
```

## References

- Hecker et al. 2022. A robust and adaptive framework for interaction
  testing in quantitative traits between multiple genetic loci and
  exposure variables.  
- Shah and Peters 2020. The Hardness of Conditional Independence Testing
  and the Generalised Covariance Measure.  
- Chernozhukov et al. 2018. Double/debiased machine learning for
  treatment and structural parameters.
- Chen. A Unified Framework for Estimation of High-dimensional
  Conditional Factor Models.
