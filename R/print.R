#' @title  print.
#' @description
#' print method for class 'jmbB'
#'
#' @param x fitted object
#' @param ... others
#'
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @importFrom rstanarm VarCorr
#' @export
#' @rdname print.jmbB
#' @method print jmbB
print.jmbB<-function(x,...){
  #x<-object
  if(!inherits(x,'jmbB'))
    stop("\n Not a 'JMbdirect' object.\n")
  result<-x
  digits<-3
  if(result$BIGdata==FALSE){
    model1<-result$model1
    model2<-result$model2
    cat('\n Joint model using jmBayes2')
    cat("\n ===========================")
    cat("\n Summary of first joint model with Event 1: ")
    cat("\n -------------------------------------------")
    cat("\n Longitudinal process:\n")
    ldat1<-list()
    ldat1$Estimate<-c(model1$statistics$Mean$betas1,model1$statistics$Mean$sigmas)
    ldat1$SE<-c(model1$statistics$SD$betas1,model1$statistics$SD$sigmas)
    ldat1$Zvalue<-ldat1$Estimate/ldat1$SE
    ldat1<-data.frame(ldat1)
    f1<-function(x){
      return(if((x[3]>0))2*(1-pnorm((x[3]),mean=0,sd=1))else 2*(pnorm((x[3]),mean=0,sd=1)))
    }
    Pvalue1<-apply(ldat1,1,f1)
    ldat1<-cbind(ldat1,Pvalue1)
    row.names(ldat1)[nrow(ldat1)]<-'sigma'
    ldat1<-round(ldat1,digits=digits)
    names(ldat1)[4]<-'Pvalue'
    print(ldat1)
    cat('\n Survival process for event1: \n ')
    sdat1<-list()
    sdat1$Estimate<-c(model1$statistics$Mean$gammas,model1$statistics$Mean$alphas)
    sdat1$SE<-c(result$model1$statistics$SD$gammas,model1$statistics$Mean$alphas)
    sdat1$Zvalue<-sdat1$Estimate/sdat1$SE
    sdat1<-data.frame(sdat1)
    Pvalue_1<-apply(sdat1,1,f1)
    sdat1<-cbind(sdat1,Pvalue_1)
    sdat1<-round(sdat1,digits=digits)
    names(sdat1)[4]<-'Pvalue'
    print(sdat1)
    cat('\n Random effect covariance matrix :\n')
    D_1<-model1$statistics$Mean$D
    print(D_1)

    cat("\n Summary of second joint model with Event 2: ")
    cat("\n -------------------------------------------")
    cat("\n Longitudinal process:\n")
    ldat2<-list()
    ldat2$Estimate<-c(model2$statistics$Mean$betas1,model2$statistics$Mean$sigmas)
    ldat2$SE<-c(model2$statistics$SD$betas1,model2$statistics$SD$sigmas)
    ldat2$Zvalue<-ldat2$Estimate/ldat2$SE
    ldat2<-data.frame(ldat2)
    f1<-function(x){
      return(if((x[3]>0))2*(1-pnorm((x[3]),mean=0,sd=1))else 2*(pnorm((x[3]),mean=0,sd=1)))
    }
    Pvalue3<-apply(ldat2,1,f1)
    ldat2<-cbind(ldat2,Pvalue3)
    row.names(ldat2)[nrow(ldat2)]<-'sigma'
    ldat2<-round(ldat2,digits=digits)
    names(ldat2)[4]<-'Pvalue'
    print(ldat2)

    cat('\n Survival process for event2: \n ')
    sdat2<-list()
    sdat2$Estimate<-c(model2$statistics$Mean$gammas,model2$statistics$Mean$alphas)
    sdat2$SE<-c(result$model2$statistics$SD$gammas,model2$statistics$Mean$alphas)
    sdat2$Zvalue<-sdat2$Estimate/sdat2$SE
    sdat2<-data.frame(sdat2)
    Pvalue_4<-apply(sdat2,1,f1)
    sdat2<-cbind(sdat2,Pvalue_4)
    sdat2<-round(sdat2,digits=digits)
    names(sdat2)[4]<-'Pvalue'
    print(sdat2)
    cat('\n Random effect covariance matrix :\n')
    D_2<-model2$statistics$Mean$D
    print(D_2)
    invisible(x)
  }else{
    model1<-result$model1
    model2<-result$model2
    cat("\n Summary of first joint model with Event 1: ")
    cat("\n -------------------------------------------")
    print(model1)
    cat("\n Summary of second joint model with Event 2: ")
    cat("\n -------------------------------------------")
    print(model2)
    invisible(x)
  }
}

#' @title  print.
#' @description
#' print method for class 'jmrmlB'
#' @param x fitted object
#' @param ... others
#'
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @importFrom rstanarm VarCorr
#' @export
#' @rdname print.jmrmlB
#' @method print jmrmlB
print.jmrmlB<-function(x,...){
  #x<-object
  if(!inherits(x,'jmrmlB'))
    stop("\n Not a 'JMbdirect' object.\n")
  result<-x
  digits<-3

  if(result$BIGdata==FALSE){
    model1<-result$model1
    model2<-result$model2

    cat("\n  Joint model using joineRML")
    cat("\n ===========================================================")
    cat("\n Summary of first joint model with Event 1: ")
    cat("\n -------------------------------------------")
    cat("\n Longitudinal process:\n")
    #cat("\n Total number of events in survival data:")

    vc1<-vcov(model1)
    nd1<-sum(1:dim(model1$coefficients$D)[[2]])
    nb1<-(nd1+1):(nd1+length(model1$coefficients$beta)+1)
    ng1<-(nd1+length(model1$coefficients$beta)+1+1):dim(vc1)[[2]]
    ese1<-diag(vc1)
    ldat1<-list()
    ldat1$Estimate<-c(model1$coefficients$beta,model1$coefficients$sigma2)
    ldat1$SE<-c(ese1[nb1])
    ldat1$Zvalue<-ldat1$Estimate/ldat1$SE
    ldat1<-data.frame(ldat1)
    f1<-function(x){
      return(if((x[3]>0))2*(1-pnorm((x[3]),mean=0,sd=1))else 2*(pnorm((x[3]),mean=0,sd=1)))
    }
    Pvalue1<-apply(ldat1,1,f1)
    ldat1<-cbind(ldat1,Pvalue=Pvalue1)
    ldat1<-round(ldat1,digits=digits)
    print(ldat1)
    cat('\n Survival process for event1: \n ')
    sdat1<-list()
    sdat1$Estimate<-model1$coefficients$gamma
    sdat1$SE<-ese1[ng1]
    sdat1$ZValue<-sdat1$Estimate/sdat1$SE
    sdat1<-data.frame(sdat1)
    sPvalue1<-apply(sdat1,1,f1)
    sdat1<-cbind(sdat1,Pvalue=sPvalue1)
    sdat1<-round(sdat1,digits=digits)
    print(sdat1)
    cat('\n Variance Covariance matrix of Random effects:\n')
    rdat1<-model1$coefficients$D
    if(dim(rdat1)[[2]]==2){
      rvar1<-all.vars(model1$formLongRandom[[1]])
      rownames(rdat1)<-c('Intercept',paste0(rvar1[[1]]))
      colnames(rdat1)<-c('Intercept',paste0(rvar1[[1]]))
    }
    rdat1<-round(rdat1,digits=digits)
    print(rdat1)

    cat("\n Summary of second joint model with Event 2:")
    cat("\n --------------------------------------------")
    cat("\n Longitudinal process:\n")
    vc2<-vcov(model2)
    nd2<-sum(1:dim(model2$coefficients$D)[[2]])
    nb2<-(nd2+1):(nd2+length(model2$coefficients$beta)+1)
    ng2<-(nd2+length(model2$coefficients$beta)+1+1):dim(vc2)[[2]]
    ese2<-diag(vc2)
    ldat2<-list()
    ldat2$Estimate<-c(model2$coefficients$beta,model2$coefficients$sigma2)
    ldat2$SE<-c(ese2[nb2])
    ldat2$Zvalue<-ldat2$Estimate/ldat2$SE
    ldat2<-as.data.frame(ldat2)
    Pvalue2<-apply(ldat2,1,f1)
    ldat2<-cbind(ldat2,Pvalue=Pvalue2)
    ldat2<-round(ldat2,digits=digits)
    print(ldat2)
    cat('\n Survival process for event2: \n ')

    sdat2<-list()
    sdat2$Estimate<-model2$coefficients$gamma
    sdat2$SE<-ese2[ng2]
    sdat2$ZValue<-sdat2$Estimate/sdat2$SE
    sdat2<-data.frame(sdat2)
    sPvalue2<-apply(sdat2,1,f1)
    sdat2<-cbind(sdat2,Pvalue=sPvalue2)
    sdat2<-round(sdat2,digits=digits)
    print(sdat2)

    cat('\n Variance Covariance matrix of Random effects:\n')
    rdat2<-model2$coefficients$D
    if(dim(rdat2)[[2]]==2){
      rvar2<-all.vars(model2$formLongRandom[[1]])
      rownames(rdat2)<-c('Intercept',paste0(rvar2[[1]]))
      colnames(rdat2)<-c('Intercept',paste0(rvar2[[1]]))
    }
    rdat2<-round(rdat2,digits=digits)
    print(rdat2)
    invisible(x)
  }else{
    model1<-result$model1
    model2<-result$model2
    cat("\n Summary of first joint model with Event 1: ")
    cat("\n -------------------------------------------")
    print(model1)
    cat("\n Summary of second joint model with Event 2: ")
    cat("\n -------------------------------------------")
    print(model2)
    invisible(x)
  }
}


#' @title  print.
#' @description
#' print method for class 'jmcsB'
#' @param x fittedobject
#' @param ... others
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @importFrom rstanarm VarCorr
#' @export
#' @rdname print.jmcsB
#' @method print jmcsB
print.jmcsB<-function(x,...){
  #x<-object
  if(!inherits(x,'jmcsB'))
    stop("\n Not a 'JMbdirect' object.\n")
  result<-x
  digits<-3

  if(x$BIGdata==FALSE){
    cat("\n  Joint model using FastJM")
    cat("\n =========================================================")
    cat("\n Summary of first joint model with Event 1: ")
    cat("\n -------------------------------------------")
    cat("\n Longitudinal process:\n")
    #cat("\n Total number of events in survival data:" ,result$model1$PropEventType[2,2] ,'\n')
    ldat1<-list()
    ldat1$Estimate<-c(result$model1_est$Par_Beta,result$model1$sigma)
    ldat1$SE<-c(result$model1$sebeta,result$model2$sesigma)
    ldat1$Zvalue<-ldat1$Estimate/ldat1$SE
    ldat1<-data.frame(ldat1)
    f1<-function(x){
      return(if((x[3]>0))2*(1-pnorm((x[3]),mean=0,sd=1))else 2*(pnorm((x[3]),mean=0,sd=1)))
    }
    Pvalue1<-apply(ldat1,1,f1)
    ldat1<-cbind(ldat1,Pvalue=Pvalue1)
    row.names(ldat1)[nrow(ldat1)]<-'sigma^2'
    ldat1<-round(ldat1,digits=digits)
    print(ldat1)
    cat('\n Survival process for event1: \n ')
    sdat1<-list()
    sdat1$Estimate<-result$model1$gamma1
    sdat1$SE<-result$model1$segamma1
    sdat1$Zvalue<-sdat1$Estimate/sdat1$SE
    sdat1<-data.frame(sdat1)
    Pvalue_1<-apply(sdat1,1,f1)
    sdat1<-cbind(sdat1,Pvalue=Pvalue_1)
    sdat1<-round(sdat1,digits=digits)
    print(sdat1)

    cat('\n Association parameters :\n')
    adat1<-list()
    adat1$Estimate<-result$model1$nu1
    adat1$SE<-result$model1$senu1
    adat1$Zvalue<-adat1$Estimate/adat1$SE
    adat1<-data.frame(adat1)
    Pvalue_a1<-apply(adat1,1,f1)
    adat1<-cbind(adat1,Pvalue=Pvalue_a1)
    adat1<-round(adat1,digits = digits)
    random1<-all.vars(result$model1$random)

    if (length(result$model1$nu1) == 2) rownames(adat1) <- c("(Intercept)_1", paste0(random[1], "_1"))

    adat1<-round(adat1,digits=digits)
    print(adat1)

    cat('\n Variance Covariance matrix of Random effects:\n')
    rdat1<-result$model1$Sig
    if(dim(rdat1)[[2]]==2){
      rvar1<-all.vars(result$model1$random)
      rownames(rdat1)<-c('Intercept',paste0(rvar1[[1]]))
      colnames(rdat1)<-c('Intercept',paste0(rvar1[[1]]))
    }
    rdat<-round(rdat1,digits=digits)
    print(rdat1)

    cat("\n Summary of second joint model with Event 2:")
    cat("\n --------------------------------------------")
    cat("\n Longitudinal process:\n")
    #cat("\n Total number of events in survival data:" ,result$model2$PropEventType[2,2],'\n')

    ldat2<-list()
    ldat2$Estimate<-c(result$model2_est$Par_Beta,result$model2$sigma)
    ldat2$SE<-c(result$model2$sebeta,result$model2$sesigma)
    ldat2$Zvalue<-ldat2$Estimate/ldat2$SE
    ldat2<-data.frame(ldat2)
    f1<-function(x){
      return(if((x[3]>0))2*(1-pnorm((x[3]),mean=0,sd=1))else 2*(pnorm((x[3]),mean=0,sd=1)))
    }
    Pvalue2<-apply(ldat2,1,f1)
    ldat2<-cbind(ldat2,Pvalue=Pvalue2)
    row.names(ldat2)[nrow(ldat2)]<-'sigma^2'
    ldat2<-round(ldat2,digits=digits)
    print(ldat2)
    cat('\n Survival process for event2: \n ')
    sdat2<-list()
    sdat2$Estimate<-result$model2$gamma1
    sdat2$SE<-result$model2$segamma1
    sdat2$Zvalue<-sdat2$Estimate/sdat2$SE
    sdat2<-data.frame(sdat2)
    Pvalue_2<-apply(sdat2,1,f1)
    sdat2<-cbind(sdat2,Pvalue=Pvalue_1)
    sdat2<-round(sdat2,digits=digits)
    print(sdat2)


    cat('\n Association parameters :\n')
    adat2<-list()
    adat2$Estimate<-result$model2$nu1
    adat2$SE<-result$model2$senu1
    adat2$Zvalue<-adat2$Estimate/adat1$SE
    adat2<-data.frame(adat2)
    Pvalue_a2<-apply(adat2,1,f1)
    adat2<-cbind(adat2,Pvalue=Pvalue_a2)
    adat2<-round(adat2,digits = digits)
    random2<-all.vars(result$model2$random)

    if (length(result$model2$nu1) == 2) rownames(adat2) <- c("(Intercept)_1", paste0(random[1], "_1"))

    adat2<-round(adat2,digits=digits)
    print(adat2)

    cat('\n Variance Covariance matrix of Random effects:\n')
    rdat2<-result$model2$Sig
    if(dim(rdat2)[[2]]==2){
      rvar2<-all.vars(result$model2$random)
      rownames(rdat2)<-c('Intercept',paste0(rvar2[[1]]))
      colnames(rdat2)<-c('Intercept',paste0(rvar2[[1]]))
    }
    rdat2<-round(rdat2,digits=digits)
    print(rdat2)
    invisible(x)
  }else{
    cat("\n Summary of second joint model with Event 1:")
    cat("\n --------------------------------------------")
    print(result$model1)
    cat("\n Summary of second joint model with Event 2:")
    cat("\n --------------------------------------------")
    print(result$model2)
    invisible(x)
  }
}


#' @title  print.
#' @description
#' print method for class 'jmstB'
#' @param x fitted object
#' @param ... others
#' @return prints table containing various parameter estimates,
#'         SE, P- value for both survival and longitudinal submodel,
#'         if the model is bayesian it includes their credible interval too.
#' @importFrom rstanarm VarCorr
#' @import dplyr
#' @export
#' @rdname print.jmstB
#' @method print jmstB
print.jmstB<-function(x,...){
  if(!inherits(x,'jmstB'))
    stop("\n Not a 'JMbdirect' object.\n")
  result<-x
  digits<-3

  model1<-result$model1
  model2<-result$model2
  if(result$BIGdata==FALSE){
  cat('\n  Joint model using rstanarm')
  cat("\n ===========================================================")
  cat("\n Summary of first joint model with Event 1: ")
  cat("\n -------------------------------------------")
  cat("\n Longitudinal process:\n")
  ldat1<-list()

  vname<-rownames(attr(model1$terms$Long1,'factors'))
  vname<-vname[vname!='id']
  vname[1]<-'Intercept'
  ldat1<-model1$stan_summary[starts_with('Long1',vars =rownames(model1$stan_summary)),c(1,3,4,10)]
  ldat1<-cbind(ldat1,Zvalue=ldat1[,1]/ldat1[,2])
  f1<-function(x){
    return(if((x[1]>0))2*(1-pnorm((x[1]/x[2]),mean=0,sd=1))else 2*(pnorm((x[1]/x[2]),mean=0,sd=1)))
  }
  attr(ldat1,'dimnames')[[2]]<-c('Mean','StDev','2.5%','97.5%','Zvalue')
  Pvalue<-apply(ldat1,1,f1)
  ldat1<-cbind(ldat1,Pvalue)
  row.names(ldat1)<-gsub("Long1\\|","",row.names(ldat1))
  #dat<-data.frame(dat,row.names = c(names(x$pseudoMod$statistics$Mean$betas1),'sigma'))
  ldat1<-round(ldat1,digits=digits)
  print(ldat1)

  cat('\n Survival process for event1: \n ')
  sdat1<-list()
  sdat1<-model1$stan_summary[c(starts_with('Event',vars =rownames(model1$stan_summary)),
                               ends_with('etavalue',vars =rownames(model1$stan_summary)))
                             ,c(1,3,4,10)]
  sdat1<-cbind(sdat1,Zvalue=sdat1[,1]/sdat1[,2])
  attr(sdat1,'dimnames')[[2]]<-c('Mean','StDev','2.5%','97.5%','Zvalue')
  sdat1<-cbind(sdat1,Pvalue=apply(sdat1,1,f1))
  row.names(sdat1)<-gsub('Event\\|','',row.names(sdat1))
  sdat1<-round(sdat1,digits=digits)
  print(sdat1)

  cat('\n Random effect covariance matrix :\n')
  D_1<-VarCorr(model1)
  print(D_1)



  cat("\n Summary of second joint model with Event 2: ")
  cat("\n -------------------------------------------")
  cat("\n Longitudinal process:\n")
  ldat2<-list()

  vname<-rownames(attr(model2$terms$Long1,'factors'))
  vname<-vname[vname!='id']
  vname[1]<-'Intercept'
  ldat2<-model2$stan_summary[starts_with('Long1',vars =rownames(model2$stan_summary)),c(1,3,4,10)]
  ldat2<-cbind(ldat2,Zvalue=ldat2[,1]/ldat2[,2])
  f1<-function(x){
    return(if((x[1]>0))2*(1-pnorm((x[1]/x[2]),mean=0,sd=1))else 2*(pnorm((x[1]/x[2]),mean=0,sd=1)))
  }
  attr(ldat2,'dimnames')[[2]]<-c('Mean','StDev','2.5%','97.5%','Zvalue')
  Pvalue<-apply(ldat2,1,f1)
  ldat1<-cbind(ldat2,Pvalue)
  row.names(ldat2)<-gsub("Long1\\|","",row.names(ldat2))
  #dat<-data.frame(dat,row.names = c(names(x$pseudoMod$statistics$Mean$betas1),'sigma'))
  ldat2<-round(ldat2,digits=digits)
  print(ldat2)

  cat('\n Survival process for event1: \n ')
  sdat2<-list()
  sdat2<-model1$stan_summary[c(starts_with('Event',vars =rownames(model2$stan_summary)),
                               ends_with('etavalue',vars =rownames(model2$stan_summary)))
                             ,c(1,3,4,10)]
  sdat2<-cbind(sdat2,Zvalue=sdat2[,1]/sdat2[,2])
  attr(sdat2,'dimnames')[[2]]<-c('Mean','StDev','2.5%','97.5%','Zvalue')
  sdat2<-cbind(sdat2,Pvalue=apply(sdat2,1,f1))
  row.names(sdat2)<-gsub('Event\\|','',row.names(sdat2))
  sdat2<-round(sdat2,digits=digits)
  print(sdat2)

  cat('\n Random effect covariance matrix :\n')
  D_2<-VarCorr(model2)
  print(D_2)
  invisible(x)
  }else{
    cat("\n Summary of first joint model with Event 1: ")
    cat("\n -------------------------------------------")
    print(model1)
    cat("\n Summary of first joint model with Event 2: ")
    cat("\n -------------------------------------------")
    print(model2)
    invisible(x)
  }
}


utils::globalVariables(c('pnorm','starts_with','ends_with','VarCorr','random','vcov'))
