#' @title Joint model for Bidirectional survival data using \code{FastJM}
#' @description
#' The function fits joint model for survival data with two events. It utilizes the FastJM package for obtaining the model parameter estimates.
#' @param dtlong longitudinal data
#' @param dtsurv survival data with two event status along with event time
#' @param longm longitudinal model e.g. list(serBilir~drug * year,serBilir ~ drug * year)
#' @param survm survival model e.g. list(Surv(years,status2)~drug,Surv(time_2,status_2)~drug+age)
#' @param rd random effect component e.g. list(~year|id,~year|id)
#' @param id ID variable
#' @param timeVar time variable
#' @param samplesize samplesize for bigdata
#' @param BIGdata logical argument TRUE or FALSE
#' @return Estimated model parameters of Joint model with bidirectional survival data
#' @importFrom FastJM jmcs
#' @import jmBIG
#' @export
#' @references
#'   Li, Shanpeng, et al. "Efficient Algorithms and Implementation of a Semiparametric Joint Model for Longitudinal and Competing Risk Data: With Applications to Massive Biobank Data." Computational and Mathematical Methods in Medicine 2022 (2022).
#'
#'   Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
#' @examples
#' library(FastJM)
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
#' pbc2_new<-pbc2[pbc2$id%in%c(1:50),]
#' new_pbc2id<-new_pbc2id[new_pbc2id$id%in%c(1:50),]
#' model_jmcs<-jmcsB(dtlong=pbc2_new,dtsurv=new_pbc2id,
#'                   longm=list(serBilir~drug*year,
#'                              serBilir~drug*year),
#'                   survm=list(Surv(years,status2)~drug,
#'                              Surv(time_2,status_2)~drug+age),
#'                   rd=list(~1|id,~1|id),
#'                   id='id',timeVar='year')
#' model_jmcs
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
jmcsB<-function(dtlong,dtsurv,longm,survm,rd,id,timeVar,BIGdata=FALSE,samplesize=200){

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
  dtlong<-as.data.frame(dtlong)
  dtsurv<-as.data.frame(dtsurv)
  if(names(dtlong)[names(dtlong)==id]=='id'){dtlong<-dtlong}else{
    dtlong<-dtlong; names(dtlong)[names(dtlong)==id]<-'id'}
  if(names(dtsurv)[names(dtsurv)==id]=='id'){dtsurv<-dtsurv}else{
    dtsurv<-dtsurv; names(dtsurv)[names(dtsurv)==id]<-'id'}
  # Preparing the data
  longm1<-longm[[1]];longm2<-longm[[2]]
  survm1<-survm[[1]];survm2<-survm[[2]]
  rd1<-rd[[1]];rd2<-rd[[2]]
  nr<-nrow(dtsurv)
  surv_st1<-all.vars(survm1)[[2]];surv_st2<-all.vars(survm2)[[2]]
  all_variable<-Reduce(union,list(all.vars(longm1),all.vars(longm2),all.vars(survm1),all.vars(survm2),all.vars(rd1),all.vars(rd2)))
  if(anyNA(dtlong[all_variable])){
    dtlong<-dtlong[complete.cases(dtlong[all_variable]),]
    dtsurv<-dtsurv[dtsurv$id%in%intersect(dtsurv$id,unique(dtlong$id)),]
  }
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
    model1<-jmcs(long.formula =longm1,ydata = data.frame(dtlong),surv.formula = survm1,
                 cdata = data.frame(dtsurv),random=rd1)
    model2<-jmcs(long.formula =longm2,ydata = data.frame(dtlong),surv.formula = survm2,
                 cdata = data.frame(dtsurv),random=rd2)
  }else{

    model1<-jmcsBig(dtlong=dtlong,dtsurv=dtsurv,longm=longm1,survm=survm1,
                    samplesize=samplesize,rd=rd1,id=id)

    model2<-jmcsBig(dtlong=dtlong,dtsurv=dtsurv,longm=longm2,survm=survm2,
                    rd=rd2,samplesize=samplesize,id=id)

  }

  result<-list()
  result$model1<-model1
  result$model1_est<-list(Par_Beta=model1$beta,Par_Gamma=model1$gamma1,Par_nu=model1$nu1,Par_vcov=model1$vcov)
  result$model2<-model2
  result$model2_est<-list(Par_Beta=model2$beta,Par_Gamma=model2$gamma1,Par_nu=model2$nu1,Par_vcov=model2$vcov)
  result$IDvar<-id
  result$timeVar<-timeVar
  result$BIGdata<-BIGdata
  class(result)<-'jmcsB'
  #samplesize=samplesize# sample of size 50 from the survival data
  #dtsurv1<-split(dtsurv,rep(1:ceiling(nr/samplesize),each=samplesize,length.out=nr))
  return(result)
}


utils::globalVariables(c('jmcs','ID','na.omit'))






