## various supporting functions
GTstat <- function(Y,A){
  ## sufficient stats
  ng = as.vector(table(A)); m = length(ng); N = sum(ng)
  sse = sum(tapply(Y, A, var)*(ng-1), na.rm=TRUE)
  mus = tapply(Y, A, mean)
  return(list(ng=ng,mus-mus,sse=sse))
}
GTpval <- function(obj, rho0=2.9){
  Qt = obj$Qt
  res = sapply(rho0^2, function(rho2){
    ans = pchisq( (rho2-Qt$Qw-Qt$Qb)/Qt$Vb,1,Qt$Mub^2/Qt$Vb, lower=FALSE )
    mean(ans)
  })
  return( res )
}
GTci <- function(obj, alpha=0.1){
  Pfc = function(x) GTpval(obj, x)-(1-alpha/2)
  x2 = 3; while(Pfc(x2)>0) x2=2*x2
  ciL = uniroot(Pfc, c(0,x2))$root
  Pfc = function(x) GTpval(obj, x)-alpha/2
  while(Pfc(x2)>0) x2=2*x2
  ciU = uniroot(Pfc, c(0,x2))$root
  return(ci=c(ciL,ciU))
}

#' Generalized test for the total absolute variation in a random-effects ANOVA model
#'
#' Given an one-way random-effects ANOVA model, \eqn{Y_{ij}=\mu+u_i+\epsilon_{ij}}
#' with \eqn{u_i\sim N(0,\sigma_b^2), \epsilon_{ij}\sim N(0,\sigma_w^2)}, we
#' want to make inference on the total absolute variation parameter, \eqn{\rho^2=\mu^2+\sigma_b^2+\sigma_w^2}.
#' A generalized inference method is taken to conduct significance test and calculate CI for \eqn{\rho}.
#'
#' @param ng sample sizes for each factor level
#' @param mus sample mean values at each factor level
#' @param sse total error sum of squares (sum of within-group sample variances)
#' @param alpha desired CI coefficient
#' @param rho0 the null threshold value to test \eqn{H_0: \rho\ge \rho_0}.
#' @param Bmc number of Monte Carlo simulations to approximate the distribution of generalized test statistics
#' @return
#' \describe{
#'   \item{p.value}{ generalized test p-value }
#'   \item{GCI}{ computed CI for \eqn{\rho} }
#'   \item{Qt}{ simulated generalized test statistics }
#' }
#'
#' @export
#' @references
#' Tsui,K.-W. and Weerahandi,S. (1989). Generalized p-Values in Significance Testing of Hypotheses in the Presence of Nuisance Parameters. Journal of the American Statistical Association 84, 602–607.
#' 
#' Weerahandi,S. (1993). Generalized Confidence Intervals. Journal of the American Statistical Association 88, 899–905.
#' 
#' Bai,Y., Wang,Z., Lystig,T., and Wu,B. (2019) Efficient and powerful equivalency test on combined mean and variance with application to diagnostic device comparison studies. arXiv:1908.07979
#' @examples
#' ## Simulation example: how to get summary stats from data
#' A = rep(1:10, 5:14)
#' Y = 0.5+rnorm(length(A))+rnorm(length(unique(A)))[A]
#' ng = as.vector(table(A))
#' sse = sum(tapply(Y, A, var)*(ng-1), na.rm=TRUE)
#' mus = tapply(Y, A, mean)
#' ans = GTrms(ng,mus,sse, alpha=0.05,rho0=2.5); ans[1:2]
#' ## Z-test
#' Zrms(ng,mus,sse, alpha=0.05,rho0=2.5)
#' ## PO comparison study
#' ng = c(9, 10, 10, 10, 5, 10, 10, 10, 10, 10, 10, 10, 2, 10, 10, 10)
#' mus = c(-0.026,0.447,0.083,-0.103,-2.587,-0.61,0.04,-0.593, 0.963,0.643,-0.2,-1.337,-4.333,-2.807,0.563,-0.797)
#' sse = 221.037
#' set.seed(123)
#' ans = GTrms(ng,mus,sse, alpha=0.1, rho0=3, Bmc=1e4); ans[1:2]
#' ## Z-test
#' Zrms(ng,mus,sse, alpha=0.1,rho0=3)
GTrms <- function(ng,mus,sse, alpha=NULL, rho0=NULL, Bmc=1e4){
  m = length(ng); N = sum(ng)
  ## inf
  Qw = sse/rchisq(Bmc,N-m); SSR = rchisq(Bmc,m-1)
  Qb = Mub = Vb = rep(0, Bmc)
  for(i in 1:Bmc){
    Qb[i] = S2Bsol(Qw[i],ng,mus,SSR[i])
    W = 1/(Qb[i] + Qw[i]/ng);  Ws = sum(W)
    Mub[i] = sum(mus*W)/Ws; Vb[i] = 1/Ws
  }
  ## GPQ
  Pfunc = function(qx){
    res = sapply(qx, function(rho2){
      ans = pchisq( (rho2-Qw-Qb)/Vb,1,Mub^2/Vb, lower=FALSE )
      mean(ans)
    })
    return(res)
  }
  ## P-value
  p.val = NULL; if(!is.null(rho0)) p.val = Pfunc(rho0^2)
  ## GCI
  ci = NULL
  if(!is.null(alpha)){
    x2 = median(Qw+Qb+Mub^2+Vb)
    fc = function(x) 1-Pfunc(x)-alpha/2
    while(fc(x2)<0) x2 = x2*2
    Lci = uniroot(fc, c(0,x2))$root
    fc = function(x) 1-Pfunc(x)-(1-alpha/2)
    while(fc(x2)<0) x2 = x2*2
    Uci = uniroot(fc, c(0,x2))$root
    ci = sqrt( c(Lci,Uci) )
  }
  ans = list(p.value=p.val, GCI=ci, Qt=list(Qw=Qw,Qb=Qb,Mub=Mub,Vb=Vb))
  class(ans) = 'RMStest'
  return( ans )
}
S2Bsol <- function(s2,ng,mus,SSR){
  fc = function(s2a){
    W = 1/(s2a+s2/ng);  yw = sum(mus*W)/sum(W)
    sum(W*(mus-yw)^2) - SSR
  }
  if(fc(0)<=0) return(0)
  x2 = 1;   while(fc(x2)>0) x2 = x2*2
  ans = uniroot(fc, c(0,x2))$root
  return(ans)
}
#' @export
print.RMStest <- function(obj){
  cat("P-value: ", obj$p.value, '\n')
  ## cat(paste0("\n  Confidence level: ", format(100*(1-alpha), digits=2), "%\n"))
  cat("Confidence interval: ", obj$GCI, '\n')
  invisible(obj)
}

#' Large-sample Z-tests for the total absolute variation in a random-effects ANOVA model
#'
#' Given a random-effects ANOVA model, \eqn{Y_{ij}=\mu+u_i+\epsilon_{ij}}
#' with \eqn{u_i\sim N(0,\sigma_b^2), \epsilon_{ij}\sim N(0,\sigma_w^2)}, we
#' want to make inference on the total absolute variation parameter, \eqn{\rho^2=\mu^2+\sigma_b^2+\sigma_w^2}.
#' Large-sample Z-tests can be conducted based on approximating \eqn{\sum_{i,j}Y_{ij}^2} with a normal distribution.
#'
#' @param ng sample sizes for each factor level
#' @param mus sample mean values at each factor level
#' @param sse total error sum of squares (sum of within-group sample variances)
#' @param alpha desired CI coefficient
#' @param rho0 the null threshold value to test \eqn{H_0: \rho\ge \rho_0}
#' @param REML whether we use REML or MLE, default to TRUE
#' @return
#' \describe{
#'   \item{p.value}{ Z-score and Z-Wald test p-values }
#'   \item{par0}{ estimated parameters, test statistic and its variance, and CI for \eqn{\rho^2}, under the null hypothesis }
#'   \item{pars}{ estimated parameters, test statistic and its variance, and CI for \eqn{\rho^2} }
#' }
#'
#' @export
#' @references
#' Ndikintum, N. K. and Rao, M. (2016). A Special Inference Problem in Repeated Measures Design Test of Statistical Hypothesis on Accuracy Root Mean Square Application to Pulse Oximetry Studies. Statistics in Biopharmaceutical Research 8, 60–76.
#'
#' Pennello, G. A. (2002). Sample size for paired repeated measures designs of method compari- son studies: Application to pulse oximetry. ASA Proceedings of the Joint Statistical Meetings pages 2671–2675.
#'
#' Pennello, G. A. (2003). Comparing monitoring devices when a gold standard in unavailable: Application to pulse oximeters. ASA Proceedings of the Joint Statistical Meetings pages 3256–3263.
#' 
#' Bai,Y., Wang,Z., Lystig,T., and Wu,B. (2019) Efficient and powerful equivalency test on combined mean and variance with application to diagnostic device comparison studies. arXiv:1908.07979
#' @examples
#' A = rep(1:10, 5:14)
#' Y = 0.5+rnorm(length(A))+rnorm(length(unique(A)))[A]
#' ng = as.vector(table(A))
#' sse = sum(tapply(Y, A, var)*(ng-1), na.rm=TRUE)
#' mus = tapply(Y, A, mean)
#' Zrms(ng,mus,sse, rho0=2.5)
Zrms <- function(ng,mus,sse, alpha=0.1, rho0=3, REML=TRUE){
  ## ng = as.vector( table(A) ); N = sum(ng); K = length(ng)
  ## mus = tapply(Y, A, mean)
  ## sse = sum( tapply(Y, A, function(yi) sum((yi-mean(yi))^2) ) )
  N = sum(ng); K = length(ng); rho2 = rho0^2
  Qy = sse + sum(ng*mus^2)
  ## score Z
  lfn = function(xpar){
    mu = xpar[1]; s2b = exp(xpar[2]); s2w = exp(xpar[3])
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*xpar[3] + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu)^2 )
    ll1 + REML*log( sum(ng/(s2w+ng*s2b)) )
  }
  cfn = function(xpar){
    mu = xpar[1]; s2b = exp(xpar[2]); s2w = exp(xpar[3])
    mu^2+s2b+s2w - rho2
  }
  options('nloptr.show.inequality.warning'=FALSE)
  xpar = nloptr::cobyla(c(mean(mus),log(rho2/2),log(rho2/2)), lfn, hin=cfn, control=list(maxeval=1e4))$par
  mu = xpar[1];  s2b = exp(xpar[2]);  s2w = exp(xpar[3])
  tau2 = ( 2*(s2w+ng*s2b)^2+2*(ng-1)*s2w^2 + 4*ng*(s2w+ng*s2b)*mu^2 )/ng^2
  V0 = sum(tau2*ng^2/N^2)
  Zs = (Qy/N-rho2)/sqrt(V0)
  pvals = pnorm(Zs)
  ci2s = Qy/N + c(1,-1)*qnorm(alpha/2)*sqrt(V0)
  ## Wald Z
  lfn = function(xpar){
    s2b = exp(xpar[1]); s2w = exp(xpar[2])
    aa1 = ng/(s2w+ng*s2b); aas = sum(aa1)
    mu = sum(aa1*mus)/aas
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*xpar[2] + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu)^2 )
    ll1 + REML*log(aas)
  }
  xpar = nloptr::newuoa(c(0,0), lfn, control=list(maxeval=1e4))$par
  s2b1 = exp(xpar[1]); s2w1 = exp(xpar[2])
  aa1 = ng/(s2w1+ng*s2b1); mu1 = sum(aa1*mus)/sum(aa1)
  tau2 = ( 2*(s2w1+ng*s2b1)^2+2*(ng-1)*s2w1^2 + 4*ng*(s2w1+ng*s2b1)*mu1^2 )/ng^2
  V1 = sum(tau2*ng^2/N^2)
  Zw = (Qy/N-rho2)/sqrt(V1)
  pvalw = pnorm(Zw)
  ci2w = Qy/N + c(1,-1)*qnorm(alpha/2)*sqrt(V1)
  ##
  pval = c(pvals,pvalw)
  names(pval) = c('Z-score', 'Z-Wald')
  return( list(p.value=pval, pars0=c(s2w=s2w,s2b=s2b,mu=mu,R=Qy/N,V=V0,ci2=ci2s), pars=c(s2w=s2w1,s2b=s2b1,mu=mu1,R=Qy/N,V=V1,ci2=ci2w)) )
}


## Random-effects ANOVA Model parameter estimation
RAMpar <- function(ng,mus,sse, mu0=0,REML=TRUE){
  N = sum(ng); K = length(ng)
  ## score Z
  lfn = function(xpar){
    s2b = exp(xpar[1]); s2w = exp(xpar[2])
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*xpar[2] + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu0)^2 )
    ll1 ## + REML*log( sum(ng/(s2w+ng*s2b)) )
  }
  options('nloptr.show.inequality.warning'=FALSE)
  xpar = nloptr::newuoa(c(0,0), lfn)$par
  s2b0 = exp(xpar[1]); s2w0 = exp(xpar[2])
  ## Wald Z
  lfn = function(xpar){
    s2b = exp(xpar[1]); s2w = exp(xpar[2])
    aa1 = ng/(s2w+ng*s2b); aas = sum(aa1)
    mu = sum(aa1*mus)/aas
    ll1 = sum( log(s2w+ng*s2b) ) + (N-K)*xpar[2] + sse/s2w + sum( ng/(s2w+ng*s2b)*(mus-mu)^2 )
    ll1 + REML*log(aas)
  }
  xpar = nloptr::newuoa(c(0,0), lfn)$par
  s2b1 = exp(xpar[1]); s2w1 = exp(xpar[2])
  aa1 = ng/(s2w1+ng*s2b1); mu1 = sum(aa1*mus)/sum(aa1)
  ##
  return( list(pars0=c(s2b=s2b0,s2w=s2w0,mu=mu0), pars=c(s2b=s2b1,s2w=s2w1,mu=mu1) ) )
}



