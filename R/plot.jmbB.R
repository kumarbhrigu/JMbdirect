#' @title Prediction plot from \code{jmbB()}
#' @param x fitted model
#' @param y newdata
#' @param ... others
#' @note
#' In the example code we use newdata as the data for ID 2 in the PBC2 dataset, it has follow up information till
#' 8.832. Now suppose we want to look at the survival of ID 2 under joint model
#' 1 after time 4 and for joint model 2 after time 9. For that we created the
#' newdata as if the individual is followed till for a time period
#' less than min(4,9).
#' @return Returns prediction plot for the newdata using the model fitted through \code{jmbB()}.
#' @import grDevices
#' @import graphics
#' @import JMbayes2
#' @import jmBIG
#' @importFrom stats predict
#' @importFrom stats median
#' @export
#'
#' @examples
#'  \donttest{
#' ##
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
#' model_jmbBdirect<-jmbB(dtlong=pbc2_new,dtsurv = new_pbc2id,
#'                        longm=list(serBilir ~ drug*year,serBilir ~ drug*year),
#'                        survm=list(Surv(years,status2) ~ drug,
#'                                   Surv(time_2,status_2) ~ drug+age),
#'                        rd=list(~year|id,~year|id),
#'                        id='id',timeVar ='year')
#'
#' t0<-4
#' nd <- pbc2[pbc2$id %in% c(2), ]
#' nd<-nd[nd$year<t0,]
#' nd$status2<-0
#' nd$years<-t0
#' nd$time_2<-9
#' nd$status_2<-0
#' plot(model_jmbBdirect,nd)
#' nd <- pbc2[pbc2$id %in% c(2), ]
#' nd<-nd[nd$year<12,]
#' nd$status2<-0
#' nd$years<-12
#' nd$time_2<-9
#' nd$status_2<-0
#' plot(model_jmbBdirect,nd)
#' ##
#'  }
#' @rdname plot.jmbB
#' @method plot jmbB
#' @export
plot.jmbB<-function(x,y,...){
  if(!inherits(x,'jmbB'))
    stop("\n Not a 'jmbBdirect' object.\n")
  result<-x
  ID<-x$IDvar
  y<-as.data.frame(y)
  newdata<-y
  t0<-max(newdata[result$timeVar])
  timeVar<-x$timeVar
  idnumber<-as.numeric(unique(newdata[,ID]))
  if(result$BIGdata==FALSE){
    long1_res<-x$model1$model_info$var_names$respVars[[1]]
    long2_res<-x$model2$model_info$var_names$respVars[[1]]
    surv1_var<-x$model1$model_info$var_names$Time_var
    surv2_var<-x$model2$model_info$var_names$Time_var
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
    pred_long2<-predJMbayes(result$model2,newdata=newdata,ids=idnumber,process='longitudinal')

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
  oldpar <- par(no.readonly = TRUE)
  K<-2
  m<-cbind(1:K,rep(K+1,K),rep(K+2,K))
  xticks<-pretty(c(pred_list[[1]]$Time,pred_list[[2]]$Time))
  widths<-c(0.3,0.35,0.35)
  layout(m,widths=widths)
  par(mar=c(0,4.5,0,0),oma=c(4,0,3,0),cex.axis=1.1,font.axis=2,font.lab=2,font.main=2)
  plot(x=pred_list[[1]]$Time,y=pred_list[[1]]$yvalue,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=pred_list[[1]]$Time,y=pred_list[[1]]$yvalue,pch=8,col='red')
  box(lwd=1.3)
  mtext(long1_res,line=2.5,side=2,font=2)
  par(mar = c(0, 4.5, 0, 0))
  plot(x=pred_list[[2]]$Time,y=pred_list[[2]]$yvalue,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=pred_list[[2]]$Time,y=pred_list[[2]]$yvalue,pch=8,col='red')
  axis(1, at = xticks)
  box(lwd=1.3)
  mtext(long2_res,line=2.5,side=2,font=2)
  xticks<-pretty(c(pred_list[[3]]$Time))
  par(mar=c(0,0,0,0))
  plot(x=pred_list[[3]]$Time,y=pred_list[[3]]$Prediction,col='black',type='l',las=1,yaxt='n',xaxs='i',xaxt='n',ylab='',ylim=c(0,1))
  polygon(c(pred_list[[3]]$Time,rev(pred_list[[3]]$Time)),c(pred_list[[3]]$UL,rev(pred_list[[3]]$LL)),col=rgb(0,0,0.6,0.4),border=NA)
  lines(x=pred_list[[3]]$Time,y=pred_list[[3]]$Prediction,col='blue',lwd=2)
  text(x=median(pred_list[[3]]$Time),y=0.99,labels =ifelse(isTRUE(survtime[1]<survtime[2]),'Event 1','Event 2'), pos = 3, col = "black",cex=1.6,font=2)
  axis(1, at = xticks)
  box(lwd=1.3)
  xticks<-pretty(c(pred_list[[4]]$Time))
  par(mar=c(0,0,0,4.5))
  plot(x=pred_list[[4]]$Time,y=pred_list[[4]]$Prediction,col='black',type='l',las=1,yaxt='n',xaxs='i',xaxt='n',ylab='',ylim=c(0,1))
  polygon(c(pred_list[[4]]$Time,rev(pred_list[[4]]$Time)),c(pred_list[[4]]$UL,rev(pred_list[[4]]$LL)),col=rgb(0.6,0,0.6,0.4),border = NA)
  lines(x=pred_list[[4]]$Time,pred_list[[4]]$Prediction,col='red',lwd=2)
  text(x=median(pred_list[[4]]$Time),y=0.99,labels = ifelse(isTRUE(survtime[1]<survtime[2]),'Event 2','Event 1'), pos = 3, col = "black",cex=1.6,font=2)
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

utils::globalVariables(c('predict'))
