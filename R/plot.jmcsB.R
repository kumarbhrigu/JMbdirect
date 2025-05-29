
#' @title Prediction plot from \code{jmcsB()}
#' @param x fitted model object
#' @param y newdata longitudinal
#' @param ... other
#' @note
#' In the example code we use newdata as the data for ID 2 in the PBC2 dataset, it has follow up information till
#' 8.832. Now suppose we want to look at the survival of ID 2 under joint model
#' 1 after time 4 and for joint model 2 after time 9. For that we created the
#' newdata as if the individual is followed till for a time period
#' less than min(4,9).
#' @return Returns prediction plot for the newdata using the model fitted through \code{jmcsB()}
#' @import jmBIG
#' @importFrom FastJM survfitjmcs
#' @importFrom stats median quantile
#' @examples
#'  \donttest{
#'  library(JMbayes2)
#'  library(FastJM)
#'   st_pbcid<-function(){
#'   new_pbcid<-pbc2.id
#'   new_pbcid$time_2<-rexp(n=nrow(pbc2.id),1/10)
#'   cen_time<-runif(nrow(pbc2.id),min(new_pbcid$time_2),max(new_pbcid$time_2))
#'   status_2<-ifelse(new_pbcid$time_2<cen_time,1,0)
#'   new_pbcid$status_2<-status_2
#'   new_pbcid$time_2<-ifelse(new_pbcid$time_2<cen_time,new_pbcid$time_2,cen_time)
#'   new_pbcid$time_2<-ifelse(new_pbcid$time_2<new_pbcid$years,new_pbcid$years,new_pbcid$time_2)
#'   new_pbcid}
#' new_pbc2id<-st_pbcid()
#' pbc2$status_2<-rep(new_pbc2id$status_2,times=data.frame(table(pbc2$id))$Freq)
#' pbc2$time_2<-rep(new_pbc2id$time_2,times=data.frame(table(pbc2$id))$Freq)
#' pbc2_new<-pbc2[pbc2$id%in%c(1:50),]
#' new_pbc2id<-new_pbc2id[new_pbc2id$id%in%c(1:50),]
#' model_jmcs<-jmcsB(dtlong=pbc2_new,dtsurv = new_pbc2id,
#'                   longm=list(serBilir~drug*year,
#'                              serBilir~drug*year),
#'                   survm=list(Surv(years,status2)~drug,
#'                              Surv(time_2,status_2)~drug+age),
#'                   rd=list(~1|id,~1|id),
#'                   id='id',timeVar='year')
#'
#' t0<-4
#' nd<-pbc2[pbc2$id %in% c(2),]
#' nd<-nd[nd$year<t0,]
#' nd$status2<-0
#' nd$years<-t0
#' nd$time_2<-9
#' nd$status_2<-0
#' plot(x=model_jmcs,y=nd)
#' ##
#'  }
#' @rdname plot.jmcsB
#' @method plot jmcsB
#' @export
plot.jmcsB<-function(x,y,...){
  if(!inherits(x,'jmcsB'))
    stop("\n Not a 'jmbBdirect' object.\n")
  result<-x
  IDvar<-x$IDvar
  y<-as.data.frame(y)
  ynewdata<-as.data.frame(y)
  cnewdata<-as.data.frame(y[!(duplicated(y[IDvar])),])
  timeVar<-x$timeVar
  idNumber<-as.numeric(cnewdata[,IDvar])
  if(result$BIGdata==FALSE){
    long1_var<-all.vars(x$model1$LongitudinalSubmodel)[1]
    long2_var<-all.vars(x$model2$LongitudinalSubmodel)[1]
    surv1_var<-all.vars(x$model1$SurvivalSubmodel)[1]
    surv2_var<-all.vars(x$model2$SurvivalSubmodel)[1]
    #future_time1<-seq(ynewdata[,timeVar][nrow(ynewdata)],3*ynewdata[,timeVar][nrow(ynewdata)],1)
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
    bootstrapped_data <- bootstrap_longitudinal_survival(x$model1$ydata,
                                                         x$model1$cdata,
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
                               long.formula=x$model1$LongitudinalSubmodel,
                               surv.formula=x$model1$SurvivalSubmodel,
                               random=x$model1$random)
        model2_list[[i]]<-jmcs(ydata=bootstrapped_data[[i]]$longitudinal,
                               cdata=bootstrapped_data[[i]]$survival,
                               long.formula=x$model2$LongitudinalSubmodel,
                               surv.formula=x$model2$SurvivalSubmodel,
                               random=x$model2$random)
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
                               long.formula=x$model2$LongitudinalSubmodel,
                               surv.formula=x$model2$SurvivalSubmodel,
                               random=x$model2$random)
        model2_list[[i]]<-jmcs(ydata=bootstrapped_data[[i]]$longitudinal,
                               cdata=bootstrapped_data[[i]]$survival,
                               long.formula=x$model1$LongitudinalSubmodel,
                               surv.formula=x$model1$SurvivalSubmodel,
                               random=x$model1$random)
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
    abline_point1<-min(pred_surv$times)
    abline_point2<-max(pred_surv$times)
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
  oldpar <- par(no.readonly = TRUE)
  K<-2
  m<-cbind(1:K,rep(K+1,K),rep(K+2,K))
  xticks<-pretty(c(y_obs$year,y_obs_2$year))
  widths<-c(0.3,0.35,0.35)
  layout(m,widths=widths)
  par(mar=c(0,4.5,0,0),oma=c(4,0,3,0),cex.axis=1.1,font.axis=2,font.lab=2,font.main=2)
  plot(x=y_obs$year,y=y_obs$Marker1,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=y_obs$year,y=y_obs$Marker1,pch=8,col='red')
  mtext(long1_var,line=2.5,side=2,font=2)
  box(lwd=1.3)
  par(mar = c(0, 4.5, 0, 0))
  plot(x=y_obs_2$year,y=y_obs_2$Marker2,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=y_obs_2$year,y=y_obs_2$Marker2,pch=8,col='red')
  mtext(long2_var,line=2.5,side=2,font=2)
  axis(1, at = xticks)
  box(lwd=1.3)
  par(mar=c(0,0,0,0))
  xticks<-pretty(c(pred_surv$times))
  plot(x=pred_surv$times,y=pred_surv$PredSurv,col='black',type='l',las=1,yaxt='n',xaxs='i',xaxt='n',ylab='',ylim=c(0,1))
  polygon(c(pred_surv$times,rev(pred_surv$times)),c(pred_surv$UL,rev(pred_surv$LL)),col=rgb(0,0,0.6,0.4),border=NA)
  lines(x=pred_surv$times,y=pred_surv$PredSurv,col='blue',lwd=2)
  text(x=median(pred_surv$times),y=0.99,labels = ifelse(t_0<=t_1,"Event 1","Event 2"), pos = 3, col = "black", cex = 1.6,font=2)
  axis(1, at = xticks)
  box(lwd=1.3)
  par(mar=c(0,0,0,4.5))
  xticks<-pretty(c(pred_surv2$times))
  plot(x=pred_surv2$times,y=pred_surv2$PredSurv,col='black',type='l',las=1,yaxt='n',xaxt='n',xaxs='i',ylab='',ylim=c(0,1))
  polygon(c(pred_surv2$times,rev(pred_surv2$times)),c(pred_surv2$UL,rev(pred_surv2$LL)),col=rgb(0.6,0,0.6,0.4),border = NA)
  lines(x=pred_surv2$times,y=pred_surv2$PredSurv,col='red',lwd=2)
  text(x=median(pred_surv2$times),y=0.99,labels = ifelse(t_0<=t_1,"Event 2","Event 1"), pos = 3, col = "black", cex = 1.6,font=2)
  axis(1, at = xticks)
  axis(4,las=1)
  mtext('Plot for Bidirectional survival data', 3,
        line = 1, outer = TRUE, font = 2, cex = 1.3)
  mtext("Time", 1,
        line = 2.5, outer = TRUE,font=2)
  mtext("Survival Prediction", 4,
        line = 2.5,font=2)
  box(lwd=1.3)
  on.exit(par(oldpar))
}

utils::globalVariables(c('survfitjmcs','year','serBilir','geom_point','times','PredSurv','geom_vline','labs','scale_y_continuous','sec_axis','scale_color_manual','theme','element_text','element_rect'))
