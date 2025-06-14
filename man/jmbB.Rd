% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jmbB.R
\name{jmbB}
\alias{jmbB}
\title{Joint model for Bidirectional survival data using \code{JMbayes2}}
\usage{
jmbB(
  dtlong,
  dtsurv,
  longm,
  survm,
  rd,
  timeVar,
  id,
  samplesize = NULL,
  BIGdata = FALSE,
  niter = 200,
  nburnin = 100,
  nchain = 1
)
}
\arguments{
\item{dtlong}{longitudinal data}

\item{dtsurv}{survival data with two event status along with event time}

\item{longm}{longitudinal model e.g. list(serBilir~drug * year,serBilir ~ drug * year)}

\item{survm}{survival model e.g. list(Surv(years,status2)~drug,Surv(time_2,status_2)~drug+age)}

\item{rd}{random effect component e.g. list(~year|id,~year|id)}

\item{timeVar}{time variable}

\item{id}{ID variable}

\item{samplesize}{samplesize for bigdata}

\item{BIGdata}{logical argument TRUE or FALSE}

\item{niter}{number if iteration}

\item{nburnin}{number of sample to burn}

\item{nchain}{number of MCMC chain}
}
\value{
Estimated model parameters of Joint model with bidirectional survival data
}
\description{
The function fits joint model for survival data with two events. It utilizes the JMbayes2 package for obtaining the model parameter estimates.
}
\examples{
library(JMbayes2)
st_pbcid<-function(){
  new_pbcid<-pbc2.id
  new_pbcid$time_2<-rexp(n=nrow(pbc2.id),1/10)
  cen_time<-runif(nrow(pbc2.id),min(new_pbcid$time_2),max(new_pbcid$time_2))
  status_2<-ifelse(new_pbcid$time_2<cen_time,1,0)
  new_pbcid$status_2<-status_2
  new_pbcid$time_2<-ifelse(new_pbcid$time_2<cen_time,new_pbcid$time_2,cen_time)
  new_pbcid$time_2<-ifelse(new_pbcid$time_2<new_pbcid$years,new_pbcid$years,new_pbcid$time_2)
  new_pbcid
}
new_pbc2id<-st_pbcid()
pbc2$status_2<-rep(new_pbc2id$status_2,times=data.frame(table(pbc2$id))$Freq)
pbc2$time_2<-rep(new_pbc2id$time_2,times=data.frame(table(pbc2$id))$Freq)
pbc2_new<-pbc2[pbc2$id\%in\%c(1:100),]
new_pbc2id<-new_pbc2id[new_pbc2id$id\%in\%c(1:100),]
model_jmbBdirect<-jmbB(dtlong=pbc2_new,dtsurv =new_pbc2id,
                       longm=list(serBilir~drug*year,serBilir~drug*year),
                       survm=list(Surv(years,status2)~drug,Surv(time_2,status_2)~drug+age),
                       rd=list(~year|id,~year|id),
                       id='id',timeVar ='year')
model_jmbBdirect
}
\references{
Rizopoulos, D., G. Papageorgiou, and P. Miranda Afonso. "JMbayes2: extended joint models for longitudinal and time-to-event data." R package version 0.2-4 (2022).

Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
}
\author{
Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
}
