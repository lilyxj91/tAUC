#' Generate a pseudo data
#'
#' @param vt Vector of time
#' @param vc Vector of censoring time
#' @param vm Vector of biomarker
#' @param vxs Matrix of biomarkers
#'
#' @return tableij
#' @import dplyr
#' @import tidyr
#' @import fastglm
#' @export
#'
data_crossing<-function(vt,vc,vm,vxs)
{
  nx<-ncol(vxs)
  nj<-length(vc)
  index<-which(vc==1)
  dat<-cbind(j=1:nj,vt,vm,vxs)
  dati<-dat[index,]%>%rename(k=j)
  nx<-ncol(vxs)
  names(dati)[1:3]<-c("k","yk","mk")
  names(dati)[4:(nx+3)]<-paste0(names(vxs),"k")
  names(dat)[1:3]<-c("j","yj","mj")
  names(dat)[4:(nx+3)]<-paste0(names(vxs),"j")
  tableij<-tidyr::crossing(dati,dat,.name_repair = "universal")
  tableij<-tableij%>%dplyr::filter(yk<yj)%>%
    dplyr::mutate(Ikj=as.numeric(mk>mj))%>%dplyr::select(-c("mk","mj"))
  tableij<-data.frame(tableij)
  return(tableij)
}
#' Standard Error estimation
#'
#' @param beta.hat coefficients
#' @param V covariance matrix
#' @param t0 time interest
#' @param nf number of elements of fractional polynomial
#'
#' @return se
#' @export
#'
se.cal<- function(beta.hat,V,t0,nf=7){
  n0 <- length(beta.hat)-nf-1
  se.store.1<- matrix(0,ncol=length(t0),nrow=3)
  for(i in 1:length(t0)){
    poly <- c(1,t0[i]^(-2),t0[i]^(-1),t0[i]^(-.5),log(t0[i]),t0[i]^(.5),t0[i],t0[i]^2,rep(0,n0))
    if (nf==3){
      poly <- c(1,t0[i]^(.5),t0[i],t0[i]^2,rep(0,n0))
    }
    if (nf==8){
      poly <- c(1,t0[i]^(-2),t0[i]^(-1),t0[i]^(-.5),log(t0[i]),t0[i]^(.5),t0[i],t0[i]^2,t0[i]^3,rep(0,n0))
    }
    temp <- c(exp(beta.hat %*% poly))
    se <- temp/(1+temp)^2*sqrt(max(0,c(poly %*% V %*% poly)) )
    se.store.1[1,i] <- se
    se.store.1[2,i] <- temp/(1+temp) - 1.96*se
    se.store.1[3,i] <- temp/(1+temp) + 1.96*se
  }
  select<-(nf+2):length(beta.hat)
  se.store.2<- matrix(0,ncol=length(select),nrow=3)
  for(i in 1:length(select)){
    se <- sqrt(max(0,diag(V)[select[i]]))
    se.store.2[1,i] <- se
    se.store.2[2,i] <- beta.hat[select[i]]- 1.96*se.store.2[1,i]
    se.store.2[3,i] <- beta.hat[select[i]] + 1.96*se.store.2[1,i]
  }
  se.store<-cbind(se.store.1,se.store.2)
  return(se.store)
}




