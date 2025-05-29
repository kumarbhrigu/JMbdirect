#' predict.jmstB
#'
#' @param object fitted model
#' @param newdata newdata
#' @param ... others
#' @return Survival Prediction for newdata from model fitted through
#' \code{jmstB()}
#' @import jmBIG
#' @export
#' @method predict jmstB
predict.jmstB<-function(object,newdata,...){
  if(!inherits(object,'jmstB'))
    stop("\n Not a 'jmbBdirect' object.\n")
  result<-object
  model1<-object$model1
  model2<-object$model2
  newdata<-as.data.frame(newdata)
  IDvar<-object$IDvar
  ID<-as.numeric(unique(newdata[,IDvar]))
  timeVar<-object$timeVar

  if(result$BIGdata==FALSE){
    first_long<-all.vars(model1$formula$Long1)[1]
    scnd_long<-all.vars(model2$formula$Long1)[1]
    id_var1<-model1$id_var
    id_var2<-model2$id_var
    long1<-model1$dataLong$Long1
    long2<-model2$dataLong$Long1
    surv1_var<-all.vars(model1$formula$Event)[1]
    surv2_var<-all.vars(model2$formula$Event)[1]
    newdata1<-newdata
    newdata2<-newdata
    t_0<-unique(newdata[,surv1_var])
    t_1<-unique(newdata[,surv2_var])
    newLdat<-newdata
    newSdat<-newdata[nrow(newdata),]
    plong1<-posterior_traj(model1,ids=ID)
    if(t_0<=t_1){
      psurv1<-posterior_survfit(model1,extrapolate=TRUE,newdataLong = newLdat,newdataEvent = newSdat,last_time =surv1_var,control=list(epoints=20) )
      psurv2<-posterior_survfit(model2,extrapolate=TRUE,newdataLong = newLdat,newdataEvent = newSdat,last_time =surv2_var,control=list(epoints=20)  )
    }else{
      psurv1<-posterior_survfit(model2,extrapolate=TRUE,newdataLong = newLdat,newdataEvent = newSdat,last_time =surv2_var,control=list(epoints=20)  )
      psurv2<-posterior_survfit(model1,extrapolate=TRUE,newdataLong = newLdat,newdataEvent = newSdat,last_time =surv1_var,control=list(epoints=20)  )

    }
    plong2<-posterior_traj(model2,ids=ID)
    if(t_0<=t_1){
      psurv1<-psurv1[psurv1$year<=t_1,]
      psurv2<-psurv2[psurv2$year<=2*t_1,]
    }else{
      psurv1<-psurv1[psurv1$year<=t_0,]
      psurv2<-psurv2[psurv2$year<=2*t_0,]
    }
    names(newdata1)[which(names(newdata1)==timeVar)]<-'year'
    names(newdata1)[which(names(newdata1)==first_long)]<-'yvalue1'
    names(newdata2)[which(names(newdata2)==timeVar)]<-'year'
    names(newdata2)[which(names(newdata2)==scnd_long)]<-'yvalue2'
  }else{
    first_long<-all.vars(model1$pseudoMod$formula$Long1)[1]
    scnd_long<-all.vars(model2$pseudoMod$formula$Long1)[1]
    id_var1<-model1$pseudoMod$id_var
    id_var2<-model2$pseudoMod$id_var
    surv1_var<-all.vars(model1$pseudoMod$formula$Event)[1]
    surv2_var<-all.vars(model2$pseudoMod$formula$Event)[1]
    newdata1<-newdata
    newdata2<-newdata
    t_0<-as.numeric(unique(newdata[,surv1_var]))
    t_1<-as.numeric(unique(newdata[,surv2_var]))
    newLdat<-newdata
    newSdat<-newdata[nrow(newdata),]
    plong1<-postTraj(model1,m=1,ids=ID)
    if(t_0<=t_1){
      psurv1<-postSurvfit(model1,ids=ID,extrapolate=TRUE,newdataLong = newLdat,
                          newdataEvent = newSdat,last_time =surv1_var,control=list(epoints=20))
      psurv2<-postSurvfit(model2,extrapolate=TRUE,ids=ID,newdataLong=newLdat,
                          newdataEvent=newSdat,last_time=surv2_var,control=list(epoints=20))
    }else{
      psurv1<-postSurvfit(model2,ids=ID,extrapolate=TRUE,newdataLong=newLdat,
                          newdataEvent=newSdat,last_time=surv2_var,control=list(epoints=20))

      psurv2<-postSurvfit(model1,ids=ID,extrapolate=TRUE,newdataLong=newLdat,
                          newdataEvent=newSdat,last_time=surv1_var,control=list(epoints=20))

    }
    plong2<-postTraj(model2,m=1,ids=ID)
    plong1<-plong1$P1
    plong2<-plong2$P1
    psurv1<-psurv1$P1
    names(psurv1)[which(names(psurv1)==timeVar)]<-'year'
    psurv2<-psurv2$P1
    names(psurv2)[which(names(psurv2)==timeVar)]<-'year'
    if(t_0<=t_1){
      psurv1<-psurv1[psurv1$year<=t_1,]
      psurv2<-psurv2[psurv2$year<=2*t_1,]
    }else{
      psurv1<-psurv1[psurv1$year<=t_0,]
      psurv2<-psurv2[psurv2$year<=2*t_0,]
    }
    names(newdata1)[which(names(newdata1)==timeVar)]<-'year'
    names(newdata1)[which(names(newdata1)==first_long)]<-'yvalue1'
    names(newdata2)[which(names(newdata2)==timeVar)]<-'year'
    names(newdata2)[which(names(newdata2)==scnd_long)]<-'yvalue2'

  }
  psurv1<-as.data.frame(psurv1)
  psurv2<-as.data.frame(psurv2)
  predlist<-list()
  predlist$long1<-plong1
  predlist$long2<-plong2
  predlist$surv1<-psurv1
  predlist$surv2<-psurv2
  predlist$ydata<-newdata1
  return(predlist)
}



#' predict.jmcsB
#'
#' @param object fitted model
#' @param newdata newdata
#' @param ... others
#' @return Survival Prediction for newdata from model fitted through
#' \code{jmcsB()}
#' @import jmBIG
#' @export
#' @method predict jmcsB
predict.jmcsB<-function(object,newdata,...){

  if(!inherits(object,'jmcsB'))
    stop("\n Not a 'jmbBdirect' object.\n")
  result<-object
  IDvar<-object$IDvar
  newdata<-as.data.frame(newdata)
  ynewdata<-as.data.frame(newdata)
  cnewdata<-as.data.frame(newdata[!(duplicated(newdata[IDvar])),])
  timeVar<-object$timeVar
  idNumber<-as.numeric(cnewdata[,IDvar])
  if(result$BIGdata==FALSE){
    long1_var<-all.vars(object$model1$LongitudinalSubmodel)[1]
    long2_var<-all.vars(object$model2$LongitudinalSubmodel)[1]
    surv1_var<-all.vars(object$model1$SurvivalSubmodel)[1]
    surv2_var<-all.vars(object$model2$SurvivalSubmodel)[1]
    t_0<-ynewdata[,surv1_var][nrow(ynewdata)]
    t_1<-ynewdata[,surv2_var][nrow(ynewdata)]
    #use bootstrap sample for 95% Confidence Interval

    bootstrap_longitudinal_survival <- function(longitudinal_data, survival_data, n_bootstrap = 10,id,idNumber){
      bootstrap_samples <- vector("list", length = n_bootstrap)
      unique_ids <- as.numeric(unique(longitudinal_data[[id]]))
      for (i in 1:n_bootstrap) {
        # Sample IDs with replacement
        bootstrap_ids <- sample(unique_ids[unique_ids!=idNumber], replace = TRUE)
        n_bootstrap_ids<-seq(1,length(unique(survival_data[[id]])))
        bootstrap_longitudinal_sample <- data.frame()
        bootstrap_survival_sample <- data.frame()
        count<-0
        for (j in 1:(length(bootstrap_ids))) {
          id_longitudinal_data <- longitudinal_data[longitudinal_data[[id]] == bootstrap_ids[j], ]
          id_longitudinal_data[[id]]<-j
          id_survival_data <- survival_data[survival_data[[id]] == bootstrap_ids[[j]], ]
          id_survival_data[[id]]<-j
          bootstrap_longitudinal_sample <- rbind(bootstrap_longitudinal_sample, id_longitudinal_data)
          bootstrap_survival_sample <- rbind(bootstrap_survival_sample, id_survival_data)
        }
        newdatalong<-longitudinal_data[longitudinal_data[id]==idNumber,]
        newdatalong[id]<-length(unique(longitudinal_data[[id]]))
        bootstrap_longitudinal_sample<-rbind(bootstrap_longitudinal_sample,newdatalong)
        newdatasurv<-survival_data[survival_data[id]==idNumber,]
        newdatasurv[id]<-length(unique(longitudinal_data[[id]]))
        bootstrap_survival_sample<-rbind(bootstrap_survival_sample,newdatasurv)
        bootstrap_samples[[i]] <- list(longitudinal = bootstrap_longitudinal_sample, survival = bootstrap_survival_sample)
      }
      return(bootstrap_samples)
    }
    bootstrapped_data <- bootstrap_longitudinal_survival(object$model1$ydata,
                                                         object$model1$cdata,
                                                         n_bootstrap = 30,id=IDvar,
                                                         idNumber = idNumber)

    if(t_0<=t_1){
      future_time1<-seq(t_0,t_1,length.out=20)
      future_time2<-seq(t_1,2*t_1,length.out=20)
      pred_model1<-survfitjmcs(object=result$model1,seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time1,method='Laplace',obs.time=timeVar)
      pred_model2<-survfitjmcs(object=result$model2,seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time2,method='Laplace',obs.time=timeVar)

      model1_list<-list();model2_list<-list()
      pred_model1_list<-list();pred_model2_list<-list();CIdata_model1<-list();
      CIdata_model2<-list()
      for(i in 1:length(bootstrapped_data)){
        model1_list[[i]]<-jmcs(ydata=bootstrapped_data[[i]]$longitudinal,
                               cdata=bootstrapped_data[[i]]$survival,
                               long.formula=object$model1$LongitudinalSubmodel,
                               surv.formula=object$model1$SurvivalSubmodel,
                               random=object$model1$random)
        model2_list[[i]]<-jmcs(ydata=bootstrapped_data[[i]]$longitudinal,
                               cdata=bootstrapped_data[[i]]$survival,
                               long.formula=object$model2$LongitudinalSubmodel,
                               surv.formula=object$model2$SurvivalSubmodel,
                               random=object$model2$random)
        pred_model1_list[[i]]<-survfitjmcs(object=model1_list[[i]],seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time1,method='Laplace',obs.time=timeVar)
        pred_model2_list[[i]]<-survfitjmcs(object=model2_list[[i]],seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time2,method='Laplace',obs.time=timeVar)
        CIdata_model1[[i]]<-pred_model1_list[[i]]$Pred[[1]]$PredSurv
        CIdata_model2[[i]]<-pred_model2_list[[i]]$Pred[[1]]$PredSurv
      }

    }else{
      future_time1<-seq(t_1,t_0,length.out=20)
      future_time2<-seq(t_0,2*t_0,length.out=20)
      pred_model1<-survfitjmcs(object=result$model2,seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time1,method='Laplace',obs.time=timeVar)
      pred_model2<-survfitjmcs(object=result$model1,seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time2,method='Laplace',obs.time=timeVar)
      model1_list<-list();model2_list<-list()
      pred_model1_list<-list();pred_model2_list<-list();CIdata_model1<-list()
      CIdata_model2<-list()
      for(i in 1:length(bootstrapped_data)){
        model1_list[[i]]<-jmcs(ydata=bootstrapped_data[[i]]$longitudinal,
                               cdata=bootstrapped_data[[i]]$survival,
                               long.formula=object$model2$LongitudinalSubmodel,
                               surv.formula=object$model2$SurvivalSubmodel,
                               random=object$model2$random)
        model2_list[[i]]<-jmcs(ydata=bootstrapped_data[[i]]$longitudinal,
                               cdata=bootstrapped_data[[i]]$survival,
                               long.formula=object$model1$LongitudinalSubmodel,
                               surv.formula=object$model1$SurvivalSubmodel,
                               random=object$model1$random)
        pred_model1_list[[i]]<-survfitjmcs(object=model2_list[[i]],seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time1,method='Laplace',obs.time=timeVar)
        pred_model2_list[[i]]<-survfitjmcs(object=model1_list[[i]],seed=100,ynewdata=ynewdata,cnewdata=cnewdata,u=future_time2,method='Laplace',obs.time=timeVar)
        CIdata_model1[[i]]<-pred_model1_list[[i]]$Pred[[1]]$PredSurv
        CIdata_model2[[i]]<-pred_model2_list[[i]]$Pred[[1]]$PredSurv

      }
    }
    CIdata_model1<-Reduce('cbind',CIdata_model1)
    CIdata_model2<-Reduce('cbind',CIdata_model2)


    y_obs <- data.frame(year = pred_model1$y.obs[[1]][,timeVar],
                        Marker1 = pred_model1$y.obs[[1]][,long1_var])

    pred_surv <- data.frame(times = pred_model1$Pred[[1]]$times,
                            PredSurv = apply(CIdata_model1,1,function(x){quantile(x,0.5)}),
                            LL=apply(CIdata_model1,1,function(x){quantile(x,0.025)}),
                            UL=apply(CIdata_model1,1,function(x){quantile(x,0.975)})
    )

    y_obs_2<-data.frame(year=pred_model2$y.obs[[1]][,timeVar],
                        Marker2=pred_model2$y.obs[[1]][,long2_var])

    pred_surv2<-data.frame(times = pred_model2$Pred[[1]]$times,
                           PredSurv =apply(CIdata_model2,1,function(x){quantile(x,0.5)}),
                           LL=apply(CIdata_model2,1,function(x){quantile(x,0.025)}),
                           UL=apply(CIdata_model2,1,function(x){quantile(x,0.975)})
    )
  }else{

    long1_var<-all.vars(result$model1$pseudoMod$LongitudinalSubmodel)[1]
    long2_var<-all.vars(result$model2$pseudoMod$LongitudinalSubmodel)[1]
    surv1_var<-all.vars(result$model1$pseudoMod$SurvivalSubmodel)[1]
    surv2_var<-all.vars(result$model2$pseudoMod$SurvivalSubmodel)[1]

    future_time1<-seq(ynewdata[,surv1_var][nrow(ynewdata)],ynewdata[,surv2_var][nrow(ynewdata)],length.out=20)
    future_time2<-seq(ynewdata[,surv2_var][nrow(ynewdata)],2*ynewdata[,surv2_var][nrow(ynewdata)],length.out=20)

    t_0<-ynewdata[,surv1_var][nrow(ynewdata)]
    t_1<-ynewdata[,surv2_var][nrow(ynewdata)]

    if(t_0<=t_1){
      future_time1<-seq(t_0,t_1,length.out=20)
      future_time2<-seq(t_1,2*t_1,length.out=20)
      pred_model1<-survfitJMCS(model=result$model1,ids=idNumber,method='Laplace',u=future_time1,obs.time = timeVar)
      pred_model2<-survfitJMCS(model=result$model2,ids=idNumber,method='Laplace',u=future_time2,obs.time=timeVar)
      bootci_model1<-bootciJMCS(pred_model1,future_time = future_time1)
      bootci_model2<-bootciJMCS(pred_model2,future_time = future_time2)
    }else{
      future_time1<-seq(t_1,t_0,length.out=20)
      future_time2<-seq(t_0,2*t_0,length.out=20)
      pred_model1<-survfitJMCS(model=result$model2,ids=ID,method='Laplace',u=future_time1,obs.time = timeVar)
      pred_model2<-survfitJMCS(model=result$model1,ids=ID,method='Laplace',u=future_time2,obs.time=timeVar)
      bootci_model1<-bootciJMCS(pred_model1,future_time = future_time1)
      bootci_model2<-bootciJMCS(pred_model2,future_time = future_time2)


    }
    y_obs <- data.frame(year = pred_model1$P1$y.obs[[1]][,timeVar],Marker1 = pred_model1$P1$y.obs[[1]][,long1_var])
    pred_surv <- data.frame(times = bootci_model1$bootCI$Times,
                            PredSurv =bootci_model1$bootCI$Med ,
                            LL=bootci_model1$bootCI$LL,
                            UL=bootci_model1$bootCI$UL)
    y_obs_2<-data.frame(year=pred_model2$P1$y.obs[[1]][,timeVar],
                        Marker2=pred_model2$P1$y.obs[[1]][,long2_var])

    pred_surv2<-data.frame(times = bootci_model2$bootCI$Times,
                           PredSurv =bootci_model2$bootCI$Med,
                           LL=bootci_model2$bootCI$LL,
                           UL=bootci_model2$bootCI$UL)

  }
  pred_surv$PredSurv<-ifelse(pred_surv$PredSurv>=1,1,pred_surv$PredSurv)
  pred_surv$LL<-ifelse(pred_surv$LL>=1,1,pred_surv$LL)
  pred_surv$UL<-ifelse(pred_surv$UL>=1,1,pred_surv$UL)
  pred_surv2$PredSurv<-ifelse(pred_surv2$PredSurv>=1,1,pred_surv2$PredSurv)
  pred_surv2$LL<-ifelse(pred_surv2$LL>=1,1,pred_surv2$LL)
  pred_surv2$UL<-ifelse(pred_surv2$UL>=1,1,pred_surv2$UL)
  predlist<-list()
  predlist$long1<-y_obs
  predlist$surv1<-pred_surv
  predlist$long2<-y_obs_2
  predlist$surv2<-pred_surv2
  predlist$ydata<-ynewdata
  return(predlist)
}


#' predict.jmrmlB
#'
#' @param object fitted model
#' @param newdata newdata
#' @param ... others
#' @return Survival Prediction for newdata from model fitted through \code{jmrmlB()}
#' @import jmBIG
#' @export
#' @method predict jmrmlB
predict.jmrmlB<-function(object,newdata,...){
  if(!inherits(object,'jmrmlB'))
    stop("\n Not a 'JMbdirect' object.\n")
  result<-object
  model1<-result$model1
  model2<-result$model2
  IDvar<-object$IDvar
  newdata<-as.data.frame(newdata)
   if(result$BIGdata==FALSE){
    longdata1<-model1$data[[1]]
    survdata1<-model1$survData
    longdata2<-model2$data[[1]]
    survdata2<-model2$survData
    surv1_var<-all.vars(object$model1$formSurv)[1]
    surv2_var<-all.vars(object$model2$formSurv)[1]
    long_var1<-all.vars(model1$formLongFixed[[1]])
    time_var1<-model1$timeVar
    long_var2<-all.vars(model2$formLongFixed[[1]])
    time_var2<-model2$timeVar
    names(longdata1)[which(names(longdata1)==model1$id)]<-'id'
    names(longdata2)[which(names(longdata2)==model2$id)]<-'id'
    names(newdata)[which(names(newdata)==IDvar)]<-'id'
    ID<-unique(newdata$id)
    ydata_1<-newdata
    ydata_2<-newdata
    t_0<-unique(newdata[,surv1_var])
    t_1<-unique(newdata[,surv2_var])
    pred_long1<-data.frame(dynLong(model1,newdata=newdata)[[1]])
    pred_long2<-data.frame(dynLong(model2,newdata=newdata)[[1]])
    if(t_0<=t_1){
      pred_surv1<-dynSurv(model1,newdata=newdata,type='simulated',u=seq(t_0,t_1,length.out=20))
      pred_surv1<-pred_surv1$pred
      pred_surv2<-dynSurv(model2,newdata=newdata,type='simulated',u=seq(t_1,2*t_1,length.out=20))
      pred_surv2<-pred_surv2$pred
    }else{
      pred_surv1<-dynSurv(model2,newdata=newdata,type='simulated',u=seq(t_1,t_0,length.out=20))
      pred_surv1<-pred_surv1$pred
      pred_surv2<-dynSurv(model1,newdata=newdata,type='simulated',u=seq(t_0,2*t_0,length.out=20))
      pred_surv2<-pred_surv2$pred
    }
    names(ydata_1)[which(names(ydata_1)==long_var1[[1]])]<-'y.value1'
    names(ydata_2)[which(names(ydata_2)==long_var2[[1]])]<-'y.value2'
    names(ydata_1)[which(names(ydata_1)==time_var1)]<-'Time'
    names(ydata_2)[which(names(ydata_2)==time_var1)]<-'Time'

  }else{
    longdata1<-list()
    for(i in 1:length(object$model1$allmodel)){
      longdata1[[i]]<-as.data.frame(object$model1$allmodel[[i]]$data)
    }
    longdata1<-Reduce('rbind',longdata1)
    survdata1<-list()
    for(i in 1:length(object$model1$allmodel)){
      survdata1[[i]]<-as.data.frame(object$model1$allmodel[[i]]$survData)
    }
    survdata1<-Reduce('rbind',survdata1)
    longdata2<-list()
    for(i in 1:length(object$model2$allmodel)){
      longdata2[[i]]<-as.data.frame(object$model2$allmodel[[i]]$data)
    }
    longdata2<-Reduce('rbind',longdata2)
    survdata2<-list()
    for(i in 1:length(object$model2$allmodel)){
      survdata2[[i]]<-as.data.frame(object$model2$allmodel[[i]]$survData)
    }
    survdata2<-Reduce('rbind',survdata2)
    surv1_var<-all.vars(object$model1$pseudoMod$formSurv)[1]
    surv2_var<-all.vars(object$model2$pseudoMod$formSurv)[1]
    long_var1<-all.vars(object$model1$pseudoMod$formLongFixed[[1]])[1]
    time_var1<-object$model1$pseudoMod$timeVar
    long_var2<-all.vars(object$model2$pseudoMod$formLongFixed[[1]])[1]
    time_var2<-object$model2$pseudoMod$timeVar
    names(longdata1)[which(names(longdata1)==model1$id)]<-'id'
    names(longdata2)[which(names(longdata2)==model2$id)]<-'id'
    names(newdata)[which(names(newdata)==IDvar)]<-'id'
    ID<-unique(newdata$id)
    ydata_1<-newdata
    ydata_2<-newdata
    t_0<-unique(newdata[,surv1_var])
    t_1<-unique(newdata[,surv2_var])
    if(t_0<=t_1){
      pred_ls1<-predJRML(model=object$model1,ids=ID,dtlong=longdata1,dtsurv = survdata1,u=seq(t_0,t_1,length.out=20))
      pred_long1<-pred_ls1$plong[[1]]$pred
      pred_surv1<-pred_ls1$psurv[[1]]$pred
      pred_ls2<-predJRML(model=object$model2,ids=ID,dtlong=longdata2,dtsurv = survdata2,u=seq(t_1,2*t_1,length.out=20))
      pred_long2<-pred_ls2$plong[[1]]$pred
      pred_surv2<-pred_ls2$psurv[[1]]$pred
    }else{
      pred_ls1<-predJRML(model=object$model2,ids=ID,dtlong=longdata2,dtsurv = survdata2,u=seq(t_1,t_0,length.out=20))
      pred_long1<-pred_ls1$plong[[1]]$pred
      pred_surv1<-pred_ls1$psurv[[1]]$pred
      pred_ls2<-predJRML(model=object$model1,ids=ID,dtlong=longdata1,dtsurv = survdata1,u=seq(t_0,2*t_0,length.out=20))
      pred_long2<-pred_ls2$plong[[1]]$pred
      pred_surv2<-pred_ls2$psurv[[1]]$pred
    }
    names(ydata_1)[which(names(ydata_1)==long_var1[[1]])]<-'y.value1'
    names(ydata_2)[which(names(ydata_2)==long_var2[[1]])]<-'y.value2'
    names(ydata_1)[which(names(ydata_1)==time_var1)]<-'Time'
    names(ydata_2)[which(names(ydata_2)==time_var1)]<-'Time'
  }

  predlist<-list()
  predlist$long1<-pred_long1
  predlist$long2<-pred_long2
  predlist$surv1<-pred_surv1
  predlist$surv2<-pred_surv2
  return(predlist)
}


#' predict.jmbB
#'
#' @param object fitted model
#' @param newdata newdata
#' @param ... others
#' @return Survival Prediction for newdata from model fitted through \code{jmbB()}
#' @import JMbayes2
#' @import jmBIG
#' @export
#' @method predict jmbB
predict.jmbB<-function(object,newdata,...){
  if(!inherits(object,'jmbB'))
    stop("\n Not a 'jmbBdirect' object.\n")
  result<-object
  ID<-object$IDvar
  newdata<-as.data.frame(newdata)
  t0<-max(newdata[result$timeVar])
  timeVar<-object$timeVar
  idnumber<-as.numeric(unique(newdata[,ID]))
  if(result$BIGdata==FALSE){
    long1_res<-object$model1$model_info$var_names$respVars[[1]]
    long2_res<-object$model2$model_info$var_names$respVars[[1]]
    surv1_var<-object$model1$model_info$var_names$Time_var
    surv2_var<-object$model2$model_info$var_names$Time_var
    surv1_time<-unique(newdata[,surv1_var])
    surv2_time<-unique(newdata[,surv2_var])
    survtime<-c(surv1_time,surv2_time)
    if(isTRUE(surv1_time<surv2_time)){
      surv1_time<-survtime[1]
      surv2_time<-survtime[2]
    }else{
      surv1_time<-survtime[2]
      surv2_time<-survtime[1]
    }
    pred_long1<-predict(result$model1,newdata=newdata,process='longitudinal')
    if(isTRUE(survtime[1]<survtime[2])){
      pred_surv1<-predict(result$model1,newdata=newdata,process='event',times=seq(surv1_time,surv2_time,length.out=20))
    }else{
      pred_surv1<-predict(result$model2,newdata=newdata,process='event',times=seq(surv1_time,surv2_time,length.out=20))
    }
    pred_long2<-predict(result$model2,newdata=newdata,process='longitudinal')
    if(isTRUE(survtime[1]<survtime[2])){
      pred_surv2<-predict(result$model2,newdata=newdata,process='event',times=seq(surv2_time,2*surv2_time,length.out=20))
    }else{
      pred_surv2<-predict(result$model1,newdata=newdata,process='event',times=seq(surv2_time,2*surv2_time,length.out=20))
    }
    pred_list<-list(long1<-data.frame(Prediction=as.vector(unlist(pred_long1$preds)),
                                      UppL=as.vector(unlist(pred_long1$upp)),
                                      LwrL=as.vector(unlist(pred_long1$low)),
                                      Time=newdata[,timeVar],
                                      yvalue=newdata[,long1_res]),
                    long2<-data.frame(Prediction=as.vector(unlist(pred_long2$preds)),
                                      UppL=as.vector(unlist(pred_long2$upp)),
                                      LwrL=as.vector(unlist(pred_long2$low)),
                                      Time=newdata[,timeVar],
                                      yvalue=newdata[,long2_res]),

                    surv1<-data.frame(Prediction=unlist(pred_surv1$pred),UL=unlist(pred_surv1$upp),LL=unlist(pred_surv1$low),Time=unlist(pred_surv1$times)),
                    surv2<-data.frame(Prediction=unlist(pred_surv2$pred),UL=unlist(pred_surv2$upp),LL=unlist(pred_surv2$low),Time=unlist(pred_surv2$times))
    )
    pred_list[[3]]$Prediction<-1-pred_list[[3]]$Prediction
    pred_list[[3]]$UL<-1-pred_list[[3]]$UL
    pred_list[[3]]$LL<-1-pred_list[[3]]$LL
    pred_list[[4]]$Prediction<-1-pred_list[[4]]$Prediction
    pred_list[[4]]$UL<-1-pred_list[[4]]$UL
    pred_list[[4]]$LL<-1-pred_list[[4]]$LL
  }else{
    long1_res<-result$model1$pseudoMod$model_info$var_names$respVars[[1]]
    long2_res<-result$model2$pseudoMod$model_info$var_names$respVars[[1]]
    surv1_var<-result$model1$pseudoMod$model_info$var_names$Time_var[[1]]
    surv2_var<-result$model2$pseudoMod$model_info$var_names$Time_var[[1]]
    surv1_time<-as.numeric(unique(newdata[,surv1_var]))
    surv2_time<-as.numeric(unique(newdata[,surv2_var]))
    survtime<-c(surv1_time,surv2_time)
    if(isTRUE(surv1_time<surv2_time)){
      surv1_time<-survtime[1]
      surv2_time<-survtime[2]
    }else{
      surv1_time<-survtime[2]
      surv2_time<-survtime[1]
    }

    pred_long1<-predJMbayes(model=result$model1,ids=idnumber,newdata=newdata,process='longitudinal')
    if(isTRUE(survtime[1]<survtime[2])){
      pred_surv1<-predJMbayes(model=result$model1,newdata=newdata,ids=idnumber,process='event',times=seq(surv1_time,surv2_time,length.out=20))
    }else{
      pred_surv1<-predJMbayes(model=result$model2,newdata=newdata,ids=idnumber,process='event',times=seq(surv1_time,surv2_time,length.out=20))
    }
    pred_long2<-predJMbayes(model=result$model2,newdata=newdata,ids=idnumber,process='longitudinal')

    if(isTRUE(survtime[1]<survtime[2])){
      pred_surv2<-predJMbayes(result$model2,newdata=newdata,ids=idnumber,process='event',times=seq(surv2_time,2*surv2_time,length.out=20))
    }else{
      pred_surv2<-predJMbayes(result$model1,newdata=newdata,ids=idnumber,process='event',times=seq(surv2_time,2*surv2_time,length.out=20))
    }
    pred_list<-list(long1<-data.frame(Prediction=as.vector(unlist(pred_long1$P1[,long1_res])),
                                      UppL=as.vector(unlist(pred_long1$P1[,paste0('upp_',long1_res)])),
                                      LwrL=as.vector(unlist(pred_long1$P1[,paste0('low_',long1_res)])),
                                      Time=as.numeric(unlist(newdata[,timeVar])),
                                      yvalue=as.numeric(unlist(newdata[,long1_res]))),
                    long2<-data.frame(Prediction=as.vector(unlist(pred_long2$P1[,long2_res])),
                                      UppL=as.vector(unlist(pred_long2$P1[,paste0('upp_',long2_res)])),
                                      LwrL=as.vector(unlist(pred_long2$P1[,paste0('low_',long2_res)])),
                                      Time=as.numeric(unlist(newdata[,timeVar])),
                                      yvalue=as.numeric(unlist(newdata[,long2_res]))),

                    surv1<-data.frame(Prediction=unlist(pred_surv1$P1[,"pred_CIF"]),UL=unlist(pred_surv1$P1[,'upp_CIF']),LL=unlist(pred_surv1$P1[,'low_CIF']),Time=unlist(pred_surv1$P1[,timeVar])),
                    surv2<-data.frame(Prediction=unlist(pred_surv2$P1[,"pred_CIF"]),UL=unlist(pred_surv2$P1[,'upp_CIF']),LL=unlist(pred_surv2$P1[,'low_CIF']),Time=unlist(pred_surv2$P1[,timeVar]))
    )
    pred_list[[3]]$Prediction<-1-pred_list[[3]]$Prediction
    pred_list[[3]]$UL<-1-pred_list[[3]]$UL
    pred_list[[3]]$LL<-1-pred_list[[3]]$LL
    pred_list[[4]]$Prediction<-1-pred_list[[4]]$Prediction
    pred_list[[4]]$UL<-1-pred_list[[4]]$UL
    pred_list[[4]]$LL<-1-pred_list[[4]]$LL
  }

  predlist<-list()
  predlist$long1<-pred_list[[1]]
  predlist$long2<-pred_list[[2]]
  predlist$surv1<-pred_list[[3]]
  predlist$surv2<-pred_list[[4]]
  predlist$ydata<-newdata
  return(predlist)
}
