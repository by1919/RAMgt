# RAMgt
 - An R package for a generalized test for a Random-effects ANOVA Model.
 - Reference
    - Bai,Y., Wang,Z., Lystig,T.C., and Wu,B. (2019) Efficient and powerful equivalency test on combined mean and variance with application to diagnostic device comparison studies. *tech rep*.
 - Sample R codes
```R
 ## install the package
 ## devtools::install_github('baolinwu/RAMgt')
library(RAMgt)
## Simulation example: how to get summary stats from data
A = rep(1:10, 5:14)
Y = 0.5+rnorm(length(A))+rnorm(length(unique(A)))[A]
ng = as.vector(table(A))
sse = sum(tapply(Y, A, var)*(ng-1), na.rm=TRUE)
mus = tapply(Y, A, mean)
ans = GTrms(ng,mus,sse, alpha=0.05,rho0=2.5); ans[1:2]
## PO comparison study
ng = c(9, 10, 10, 10, 5, 10, 10, 10, 10, 10, 10, 10, 2, 10, 10, 10)
mus = c(-0.026,0.447,0.083,-0.103,-2.587,-0.61,0.04,-0.593, 0.963,0.643,-0.2,-1.337,-4.333,-2.807,0.563,-0.797)
sse = 221.037
set.seed(123)
GTrms(ng,mus,sse, alpha=0.1, rho0=3, Bmc=1e4)
## Z-test
Zrms(ng,mus,sse, alpha=0.1,rho0=3)
```
