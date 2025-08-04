
<!-- README.md is generated from README.Rmd. Please edit that file -->



# ROMY: Robust Omics MethodologY

<!-- badges: start -->
<!-- badges: end -->

The *romy* R package provides flexible and robust tools for the analysis of multi-omics data. **!This package is still under development and not for use yet.!**

## Installation

You can install *romy* using:

``` {.r}
#library(devtools)
#install_github("julianhecker/romy")
library(romy)
```

## Background

The *romy* R package implements five different analysis components for the analysis of multi-omics data.
- Conditional Independence Testing (CIT): testing for association between different omics factors and phenotype variables while adjusting for covariates/confounders Z.
- (co)variance testing (COV): testing for an effect of X on the covariance between two factors in Y, or the variance of one factor in Y, adjusting for Z.
- interaction testing (INTER): testing for interaction between factor in X and factor in Z in their effect on Y.
- Factor Conditional Independence Testing (FCIT): testing for association between an latent factor of Y and X, adjusting for Z.
- Latent Factor (LF): Estimation of latent factors that correspond to heterogenous effect profiles.

The first four approaches utilize the concepts of debiased machine learning and/or robust statistics (Chernozhukov et al. 2018, Shah and Peters 2020).
The last component (LF) is based on the work of Chen (2022).

## CIT

The first analysis approach in *romy* performs association testing between variables in Y and X, adjusting for Z. The measurements in X and Y can be transformed using data transformation functions. In this case, all pairs of transformations are tested separately and jointly. The separate p-values are also combined in an ACAT (Aggregated Cauchy Association Test).

``` {.r}
X=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Z=matrix(rnorm(1000*3), nrow=1000, ncol=3)
Y=matrix(rnorm(1000*30), nrow=1000, ncol=30)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

at=romy_cit(Y, X, Z)
head(at)
#>   X_id Y_id   p_value
#> 1   X1   Y1 0.9672532
#> 2   X2   Y1 0.9622144
#> 3   X3   Y1 0.6385682
#> 4   X1   Y2 0.4666252
#> 5   X2   Y2 0.2802362
#> 6   X3   Y2 0.4462549
```

## COV

The second analysis approach aims to test for effects of variables in X on the (co)variance in Y, adjusting for Z.

``` {.r}
X=matrix(rnorm(1000*3),nrow=1000,ncol=3)
Z=matrix(rnorm(1000*3),nrow=1000,ncol=3)
Y=matrix(rnorm(1000*10),nrow=1000,ncol=10)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

ct=romy_cov(Y, X, Z)
head(ct)
#>   X_id Y_id1 Y_id2    p_value
#> 1   X1    Y1    Y1 0.07553816
#> 2   X2    Y1    Y1 0.80649030
#> 3   X3    Y1    Y1 0.78876153
#> 4   X1    Y1    Y2 0.34012227
#> 5   X2    Y1    Y2 0.52967296
#> 6   X3    Y1    Y2 0.75683554
```

## INTER

The third analysis approach aims to perform robust interaction testing between variables in X and Z in their effect on Y. This is similar to the RITSS (Robust Interaction Testing using Sample Splitting) method.

``` {.r}
X=matrix(rnorm(1000*30),nrow=1000,ncol=30)
Z=matrix(rnorm(1000*5),nrow=1000,ncol=5)
Y=matrix(rnorm(1000*1),nrow=1000,ncol=1)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))
colnames(Z)=paste0("Z",1:ncol(Z))

it=romy_inter(Y, X, Z)
head(it)
#>   X_id Y_id   p_value
#> 1   X1   Z1 0.8209516
#> 2   X2   Z1 0.6272416
#> 3   X3   Z1 0.8296387
#> 4   X4   Z1 0.3573190
#> 5   X5   Z1 0.9853387
#> 6   X6   Z1 0.7305682
```

## FCIT

The fourth analysis approach aims to perform association analysis between X and a latent factor underlying Y, adjusting for Z.

``` {.r}
X=matrix(rnorm(1000*1),nrow=1000,ncol=1)
Z=matrix(rnorm(1000*2),nrow=1000,ncol=2)
Y=matrix(rnorm(1000*10),nrow=1000,ncol=10)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

fcit=romy_fcit(Y=Y, X=X, Z=Z)
head(fcit)
#>   X_id   p_value
#> 1   X1 0.8462599
```

## LF

The fifth analysis approach aims to identify and estimate latent factors that correspond to effect profile heterogeneity in the data.

``` {.r}
N=1000
M=30
P_m=2
P= P_m*M
K=1
beta=2.0
X=matrix(rnorm(N*P),nrow=N,ncol=P)
Y=matrix(rnorm(N*M),nrow=N,ncol=M)
colnames(X)=paste0("X",1:ncol(X))
colnames(Y)=paste0("Y",1:ncol(Y))

inds_X=list()
for(m in 1:M)
{
    inds_X[[m]]=c(((m-1)*P_m+1):((m-1)*P_m+P_m))
}
Z1=matrix(0, nrow=N, ncol=K+1)
F1=matrix(0, nrow=K+1, ncol=P)
        
for(m in 1:M)
{
          alpha_m=rnorm(P_m, mean=0, sd=1)
          F1[1,((m-1)*P_m+1):((m-1)*P_m+P_m)]=alpha_m
          for(k in 1:K)
          {
            beta_km=rnorm(P_m, mean=0, sd=1)
            F1[k+1,((m-1)*P_m+1):((m-1)*P_m+P_m)]=beta_km
          }
}
        
Z1[,1]=1
Z=matrix(rnorm(N*K, mean = 0, sd=1), nrow=N, ncol=K)
Z1[,2:(K+1)]=Z
F=F1[2:(K+1),]
        
        
Gamma=Z1 %*% F1
        
        
        
for(m in 1:M)
{
          Y[,m]=Y[,m]+beta*rowSums(as.matrix(X[,inds_X[[m]]]) * as.matrix(Gamma[,inds_X[[m]]]))
}
lf=romy_lf(Y=Y, X=X, inds_X=inds_X)
cor(Z[,1], lf$Z_hat[,1])
#> [1] -0.9821306
```

## References

-   Hecker et al. 2022. A robust and adaptive framework for interaction testing in quantitative traits between multiple genetic loci and exposure variables.
-   Shah and Peters 2020. The Hardness of Conditional Independence Testing and the Generalised Covariance Measure.
-   Chernozhukov et al. 2018. Double/debiased machine learning for treatment and structural parameters.
-   Chen. A Unified Framework for Estimation of High-dimensional Conditional Factor Models.
