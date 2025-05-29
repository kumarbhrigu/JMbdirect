
#' @title Joint model for Bidirectional survival data using \code{JMbayes2}
#' @description
#' The function fits joint model for survival data with two events. It utilizes the JMbayes2 package for obtaining the model parameter estimates.
#' @param dtlong longitudinal data
#' @param dtsurv survival data with two event status along with event time
#' @param longm longitudinal model e.g. list(serBilir~drug * year,serBilir ~ drug * year)
#' @param survm survival model e.g. list(Surv(years,status2)~drug,Surv(time_2,status_2)~drug+age)
#' @param rd random effect component e.g. list(~year|id,~year|id)
#' @param timeVar time variable
#' @param id ID variable
#' @param niter number if iteration
#' @param nburnin number of sample to burn
#' @param nchain number of MCMC chain
#' @param samplesize samplesize for bigdata
#' @param BIGdata logical argument TRUE or FALSE
#' @return Estimated model parameters of Joint model with bidirectional survival data
#' @importFrom JMbayes2 jm
#' @importFrom stats predict complete.cases
#' @import jmBIG
#' @export
#' @references
#' Rizopoulos, D., G. Papageorgiou, and P. Miranda Afonso. "JMbayes2: extended joint models for longitudinal and time-to-event data." R package version 0.2-4 (2022).
#'
#' Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
#' @examples
#' library(JMbayes2)
#' st_pbcid<-function(){
#'   new_pbcid<-pbc2.id
#'   new_pbcid$time_2<-rexp(n=nrow(pbc2.id),1/10)
#'   cen_time<-runif(nrow(pbc2.id),min(new_pbcid$time_2),max(new_pbcid$time_2))
#'   status_2<-ifelse(new_pbcid$time_2<cen_time,1,0)
#'   new_pbcid$status_2<-status_2
#'   new_pbcid$time_2<-ifelse(new_pbcid$time_2<cen_time,new_pbcid$time_2,cen_time)
#'   new_pbcid$time_2<-ifelse(new_pbcid$time_2<new_pbcid$years,new_pbcid$years,new_pbcid$time_2)
#'   new_pbcid
#' }
#' new_pbc2id<-st_pbcid()
#' pbc2$status_2<-rep(new_pbc2id$status_2,times=data.frame(table(pbc2$id))$Freq)
#' pbc2$time_2<-rep(new_pbc2id$time_2,times=data.frame(table(pbc2$id))$Freq)
#' pbc2_new<-pbc2[pbc2$id%in%c(1:100),]
#' new_pbc2id<-new_pbc2id[new_pbc2id$id%in%c(1:100),]
#' model_jmbBdirect<-jmbB(dtlong=pbc2_new,dtsurv =new_pbc2id,
#'                        longm=list(serBilir~drug*year,serBilir~drug*year),
#'                        survm=list(Surv(years,status2)~drug,Surv(time_2,status_2)~drug+age),
#'                        rd=list(~year|id,~year|id),
#'                        id='id',timeVar ='year')
#' model_jmbBdirect
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
jmbB<-function(dtlong,dtsurv,longm,survm,rd,timeVar,id,
                 samplesize=NULL,BIGdata=FALSE,niter=200,nburnin=100,nchain=1){
  cl<-match.call()
  if(!id%in%names(dtlong) ){
    stop("\n Longitudinal data must have column 'id' ")
  }
  if(!id%in%names(dtsurv) ){
    stop("\n Survival data must have column 'id' ")
  }
  if(!names(dtlong)[names(dtlong)==id]==names(dtsurv)[names(dtsurv)==id]){
    stop("\n'dtlong' and 'dtsurv' must have same id.")
  }
  if(!timeVar%in%names(dtlong)){
    stop("\n 'timeVar' should be in longitudinal dataset")
  }

  dtlong<-as.data.frame(dtlong)
  dtsurv<-as.data.frame(dtsurv)
  if(names(dtlong)[names(dtlong)==id]=='id'){dtlong<-dtlong}else{
    dtlong<-dtlong; names(dtlong)[names(dtlong)==id]<-'id'}
  if(names(dtsurv)[names(dtsurv)==id]=='id'){dtsurv<-dtsurv}else{
    dtsurv<-dtsurv; names(dtsurv)[names(dtsurv)==id]<-'id'}

  longm1<-longm[[1]];longm2<-longm[[2]]
  survm1<-survm[[1]];survm2<-survm[[2]]
  rd1<-rd[[1]];rd2<-rd[[2]]
  all_variable<-Reduce(union,list(all.vars(longm1),all.vars(longm2),all.vars(survm1),all.vars(survm2),all.vars(rd1),all.vars(rd2)))
  if(anyNA(dtlong[all_variable])){
    dtlong<-dtlong[complete.cases(dtlong[all_variable]),]
    dtsurv<-dtsurv[dtsurv$id%in%intersect(dtsurv$id,unique(dtlong$id)),]
  }
  surv_st1<-all.vars(survm1)[[2]];surv_st2<-all.vars(survm2)[[2]]
  if(is.factor(dtsurv[,surv_st1])){
    if(!is.numeric(levels(dtsurv[,surv_st1]))){
      stop('Use status variable a numeric with censored=0 and dead=1')
    }else{

      if(length(levels(dtsurv[,surv_st1]))!=2){
        stop('More than 2 possible values for survival status. Use status variable a numeric with censored=0 and dead=1')
      }

      if(length(levels(dtsurv[,surv_st1]))==2&sum(levels(dtsurv[,surv_st1]))!=1){
        stop('More than two possible survival status.Use status variable a numeric with censored=0 and dead=1')
      }
    }
  }

  if(is.factor(dtsurv[,surv_st2])){
    if(!is.numeric(levels(dtsurv[,surv_st2]))){
      stop('Use status variable as numeric with censored=0 and dead=1')
    }else{

      if(length(levels(dtsurv[,surv_st2]))!=2){
        stop('More than 2 possible values for survival status. Use status variable a numeric with censored=0 and dead=1')
      }

      if(length(levels(dtsurv[,surv_st2]))==2&sum(levels(dtsurv[,surv_st2]))!=1){
        stop('Use status variable as numeric with censored=0 and dead=1')
      }
    }
  }

  if(is.numeric(dtsurv[,surv_st1])&&length(unique(dtsurv[,surv_st1]))!=2){
    stop('More than two possible survival status.')
  }

  if(BIGdata==FALSE){
    fm1<-lme(fixed=longm1,random=rd1,data=dtlong)
    fm2<-coxph(survm1,data=dtsurv)
    mod1 <- jm(fm2,fm1,time_var=timeVar,n_chains=nchain,n_iter=niter,n_burnin=nburnin)
    fm3<-lme(fixed=longm2,random=rd2,data=dtlong)
    fm4<-coxph(survm2,data=dtsurv)
    mod2 <- jm(fm4,fm3,time_var=timeVar,n_chains=nchain,n_iter=niter,n_burnin=nburnin)
  }else{
    mod1<-jmbayesBig(dtlong=dtlong,dtsurv =dtsurv,longm=longm1,survm=survm1,
                     rd=rd1,timeVar=timeVar,nchain=nchain,samplesize=samplesize,id=id)

    mod2<-jmbayesBig(dtlong=dtlong,dtsurv =dtsurv,longm=longm2,survm=survm2,rd=rd2,
                     timeVar=timeVar,nchain=nchain,samplesize=samplesize,id=id)

  }
  result<-list()
  result$call<-cl
  result$model1<-mod1
  result$model2<-mod2
  result$timeVar<-timeVar
  result$IDvar<-id
  result$BIGdata<-BIGdata
  class(result)<-'jmbB'
  #message('Results for the joint model')
  result
}
utils::globalVariables(c('lme','coxph','jm','na.omit'))
