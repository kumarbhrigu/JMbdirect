

#' @title Joint model for Bidirectional survival data using \code{joineRML}
#' @description
#' The function fits joint model for survival data with two events. It utilizes the joineRML package for obtaining the model parameter estimates.
#' @param dtlong longitudinal data
#' @param dtsurv survival data with two event status along with event time
#' @param longm longitudinal model e.g. list(serBilir~drug * year,serBilir ~ drug * year)
#' @param survm survival model e.g. list(Surv(years,status2)~drug,Surv(time_2,status_2)~drug+age)
#' @param rd random effect component e.g. list(~year|id,~year|id)
#' @param timeVar time variable
#' @param id ID variable
#' @param samplesize samplesize for bigdata
#' @param BIGdata logical argument TRUE or FALSE
#' @return Estimated model parameters of Joint model with bidirectional survival data
#' @import joineRML
#' @export
#' @references
#' Hickey, Graeme L., et al. "joineRML: a joint model and software package for time-to-event and multivariate longitudinal outcomes." BMC medical research methodology 18 (2018): 1-14.
#'
#' Bhattacharjee, A., Rajbongshi, B. K., & Vishwakarma, G. K. (2024). jmBIG: enhancing dynamic risk prediction and personalized medicine through joint modeling of longitudinal and survival data in big routinely collected data. BMC Medical Research Methodology, 24(1), 172.
#' @examples
#'  \donttest{
#' ##
#' library(JMbayes2)
#' library(joineRML)
#' jmrmlBModel<-jmrmlB(dtlong=new_long2[new_long2$id%in%c(1:80),],
#'                     dtsurv=new_surv2[new_surv2$id%in%c(1:80),],
#'                     longm=list(y~x7+visit,y~x7+visit),survm=list(Surv(time,status)~x1+visit,
#'                     Surv(time_2,status_2)~x1+visit),rd=list(~visit|id,~visit|id),id='id',
#'                     timeVar='visit',samplesize=40,BIGdata=TRUE)
#' jmrmlBModel
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
jmrmlB<-function(dtlong,dtsurv,longm,survm,rd,timeVar,
                 id,samplesize=200,BIGdata=FALSE){
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
    model1<-mjoint(formLongFixed=longm1,formLongRandom=rd1,formSurv=survm1,data=dtlong
                   ,survData=dtsurv,timeVar= timeVar)
    model2<-mjoint(formLongFixed=longm2,formLongRandom=rd2,formSurv=survm2,data=dtlong
                   ,survData=dtsurv,timeVar= timeVar)
  }else{
    model1<-joinRMLBig(dtlong=dtlong,dtsurv = dtsurv,longm=longm1,survm=survm1,
                       rd=rd1,timeVar=timeVar,samplesize=samplesize,id=id)

    model2<-joinRMLBig(dtlong=dtlong,dtsurv = dtsurv,longm=longm2,survm=survm2,
                       rd=rd2,timeVar=timeVar,samplesize=samplesize,id=id)
  }
  result<-list()
  result$model1<-model1
  result$model2<-model2
  result$timeVar<-timeVar
  result$IDvar<-id
  result$BIGdata<-BIGdata
  class(result)<-'jmrmlB'
  result
}

utils::globalVariables(c('mjoint','na.omit'))










