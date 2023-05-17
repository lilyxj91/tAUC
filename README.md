
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tAUC

<!-- badges: start -->
<!-- badges: end -->

The goal of tAUC is to estimate the covariate-specific incident dynamic
time-dependent AUC.

## Installation

You can install the development version of tAUC like so:

``` r
library(devtools)
devtools::install_github("lilyxj91/tAUC",force = TRUE)
```

## Example

``` r
library(tAUC)
```

``` r

set.seed(1)
t0 = c(1.8, 3.6, 6.9) #day 1
head(simdata)
#>            Y delta xd         xc          M
#> 1 0.05592796     0  1 0.90729142 -2.8112120
#> 2 0.12588292     1  0 0.99278809  0.3783830
#> 3 0.12999906     0  1 0.04965358  0.6373297
#> 4 0.14487118     0  0 0.82027999 -1.5550203
#> 5 0.14498653     1  0 0.96467367  0.8961479
#> 6 0.16240420     1  0 0.40827776  0.3466833
```

``` r
Y0<-simdata$Y
C0<-simdata$delta
M0<-simdata$M
VXS<-simdata[c("xd","xc")]
datij<-data_crossing(Y0,C0,M0,VXS)
print(dim(datij))
#> [1] 99965     9
YK = datij$yk
XS<-cbind(int = 1,
          t1 = YK^(-2),
          t2 = YK^(-1),
          t3 = YK^(-0.5),
          t4 = log(YK),
          t5 = YK^(0.5),
          t6 = YK,
          t7 = YK^(2),
          xd10 = ifelse(datij$xdk==1&datij$xdj==0,1,0),
          xd01 = ifelse(datij$xdk==0&datij$xdj==1,1,0),
          xd11 = ifelse(datij$xdk==1&datij$xdj==1,1,0),
          xci = datij$xck,
          xcj = datij$xcj)
YS<-datij$Ikj
m<-fastglm::fastglm(XS,YS,family=binomial(link="logit"),maxit=10000L)
beta.hat<-m$coefficients
beta.hat.info<-c(m$converged,m$iter)
```

``` r
ordic=datij$k-1
ordjc=datij$j-1
L=cov_cal_long(beta.hat,M0,YS,t(XS),ordic,ordjc) #1 min
V=solve(L$sigma1)%*%(L$sigma2)%*%solve(L$sigma1)
cov=se.cal(beta.hat,V,t0,nf=7)
```

``` r
select=c(9:length(beta.hat))
se.store = diag(V)[select]
for(i in 1:length(select)){
  se <- sqrt(max(0.00000001,diag(V)[select[i]]))
  se.store[i] <- se
}
beta = beta.hat[select]
wald = beta/se.store
pvalue = 2*(1-pnorm(abs(wald)))
table<-round(as.matrix(cbind(
  beta,
  se.store,
  wald,
  pvalue)),3)
colnames(table)<-c("Estimate","SE","Wald","P-value")
table[,4]<-ifelse(table[,4]<0.001,"<0.001",table[,4])
table
#>      Estimate SE      Wald     P-value 
#> xd10 "-0.361" "0.144" "-2.501" "0.012" 
#> xd01 "0.983"  "0.15"  "6.557"  "<0.001"
#> xd11 "0.457"  "0.129" "3.536"  "<0.001"
#> xci  "-0.067" "0.203" "-0.329" "0.742" 
#> xcj  "0.043"  "0.237" "0.18"   "0.857"
```
