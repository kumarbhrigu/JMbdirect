#' @title Prediction plot from \code{jmrmlB()}
#' @param x fitted model object
#' @param y newdata
#' @param ... others
#' @note
#' In the example code we use newdata as the data for ID 2 in the PBC2 dataset, it has follow up information till
#' 8.832. Now suppose we want to look at the survival of ID 2 under joint model
#' 1 after time 4 and for joint model 2 after time 9. For that we created the
#' newdata as if the individual is followed till for a time period
#' less than min(4,9).
#' @return Returns prediction plot for the newdata using the model fitted through \code{jmrmlB()}.
#' @import joineRML
#' @import grDevices
#' @import graphics
#' @importFrom stats median
#' @export
#' @examples
#'  \donttest{
#' ##
#' library(JMbayes2)
#' library(joineRML)
#' jmrmlBModel<-jmrmlB(dtlong=new_long2[new_long2$id<=400,],
#'                     dtsurv=new_surv2[new_surv2$id<=400,],
#'                     longm=list(y~x7+visit,y~x7+visit),
#'                     survm=list(Surv(time,status)~x1+visit,
#'                                Surv(time_2,status_2)~x1+visit),
#'                     rd=list(~visit|id,~visit|id),
#'                     id='id',
#'                     timeVar='visit',
#'                     samplesize=200,
#'                     BIGdata=TRUE)
#' t0<-6
#' ndBIG<-new_long2[new_long2$id==10,]
#' ndBIG<-ndBIG[ndBIG$visit<t0,]
#' ndBIG$status<-0
#' ndBIG$time<-t0
#' ndBIG$time_2<-10
#' ndBIG$status_2<-0
#' plot(jmrmlBModel,ndBIG)
#' ##
#' }
#' @rdname plot.jmrmlB
#' @method plot jmrmlB
#' @export
plot.jmrmlB<-function(x,y,...){
  if(!inherits(x,'jmrmlB'))
    stop("\n Not a 'JMbdirect' object.\n")
  result<-x
  model1<-result$model1
  model2<-result$model2
  IDvar<-x$IDvar
  y<-as.data.frame(y)
  newdata<-as.data.frame(y)
  if(result$BIGdata==F){
    longdata1<-model1$data[[1]]
    survdata1<-model1$survData
    longdata2<-model2$data[[1]]
    survdata2<-model2$survData
    surv1_var<-all.vars(x$model1$formSurv)[1]
    surv2_var<-all.vars(x$model2$formSurv)[1]
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
    for(i in 1:length(x$model1$allmodel)){
      longdata1[[i]]<-as.data.frame(x$model1$allmodel[[i]]$data)
    }
    longdata1<-Reduce('rbind',longdata1)
    survdata1<-list()
    for(i in 1:length(x$model1$allmodel)){
      survdata1[[i]]<-as.data.frame(x$model1$allmodel[[i]]$survData)
    }
    survdata1<-Reduce('rbind',survdata1)
    longdata2<-list()
    for(i in 1:length(x$model2$allmodel)){
      longdata2[[i]]<-as.data.frame(x$model2$allmodel[[i]]$data)
    }
    longdata2<-Reduce('rbind',longdata2)
    survdata2<-list()
    for(i in 1:length(x$model2$allmodel)){
      survdata2[[i]]<-as.data.frame(x$model2$allmodel[[i]]$survData)
    }
    survdata2<-Reduce('rbind',survdata2)
    surv1_var<-all.vars(x$model1$pseudoMod$formSurv)[1]
    surv2_var<-all.vars(x$model2$pseudoMod$formSurv)[1]
    long_var1<-all.vars(x$model1$pseudoMod$formLongFixed[[1]])[1]
    time_var1<-x$model1$pseudoMod$timeVar
    long_var2<-all.vars(x$model2$pseudoMod$formLongFixed[[1]])[1]
    time_var2<-x$model2$pseudoMod$timeVar
    names(longdata1)[which(names(longdata1)==model1$id)]<-'id'
    names(longdata2)[which(names(longdata2)==model2$id)]<-'id'
    names(newdata)[which(names(newdata)==IDvar)]<-'id'
    ID<-unique(newdata$id)
    ydata_1<-newdata
    ydata_2<-newdata
    t_0<-unique(newdata[,surv1_var])
    t_1<-unique(newdata[,surv2_var])
    if(t_0<=t_1){
      pred_ls1<-predJRML(model=x$model1,ids=ID,dtlong=longdata1,dtsurv = survdata1,u=seq(t_0,t_1,length.out=20))
      pred_long1<-pred_ls1$plong[[1]]$pred
      pred_surv1<-pred_ls1$psurv[[1]]$pred
      pred_ls2<-predJRML(model=x$model2,ids=ID,dtlong=longdata2,dtsurv = survdata2,u=seq(t_1,2*t_1,length.out=20))
      pred_long2<-pred_ls2$plong[[1]]$pred
      pred_surv2<-pred_ls2$psurv[[1]]$pred
    }else{
      pred_ls1<-predJRML(model=x$model2,ids=ID,dtlong=longdata2,dtsurv = survdata2,u=seq(t_1,t_0,length.out=20))
      pred_long1<-pred_ls1$plong[[1]]$pred
      pred_surv1<-pred_ls1$psurv[[1]]$pred
      pred_ls2<-predJRML(model=x$model1,ids=ID,dtlong=longdata1,dtsurv = survdata1,u=seq(t_0,2*t_0,length.out=20))
      pred_long2<-pred_ls2$plong[[1]]$pred
      pred_surv2<-pred_ls2$psurv[[1]]$pred
    }
    names(ydata_1)[which(names(ydata_1)==long_var1[[1]])]<-'y.value1'
    names(ydata_2)[which(names(ydata_2)==long_var2[[1]])]<-'y.value2'
    names(ydata_1)[which(names(ydata_1)==time_var1)]<-'Time'
    names(ydata_2)[which(names(ydata_2)==time_var1)]<-'Time'
  }
  oldpar <- par(no.readonly = TRUE)
  K<-2
  m<-cbind(1:K,rep(1+K,K),rep(K+2,K))
  widths<-c(0.3,0.35,0.35)
  layout(m,widths=widths)
  xticks<-pretty(c(ydata_1$Time,ydata_2$Time))
  par(mar=c(0,4.5,0,0),oma=c(4,0,3,0),cex.axis=1.1,font.axis=2,font.lab=2,font.main=2)
  plot(x=ydata_1$Time,y=ydata_1$y.value1,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=ydata_1$Time,y=ydata_1$y.value1,pch=8,col='red')
  mtext(long_var1[1],line=2.5,side=2,font=2)
  box(lwd=1.3)
  par(mar=c(0,4.5,0,0))
  plot(x=ydata_2$Time,y=ydata_2$y.value2,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=ydata_2$Time,y=ydata_2$y.value2,pch=8,col='red')
  mtext(long_var2[1],line=2.5,side=2,font=2)
  axis(1, at = xticks)
  box(lwd=1.3)
  xticks<-pretty(c(pred_surv1$u))
  par(mar=c(0,0,0,0))
  plot(x=pred_surv1$u,y=pred_surv1$mean,col='black',type='l',las=1,yaxt='n',xaxt='n',xaxs='i',ylab='',ylim=c(0,1),font=2)
  polygon(c(pred_surv1$u,rev(pred_surv1$u)),c(pred_surv1$upper,rev(pred_surv1$lower)),col=rgb(0,0,0.6,0.4),border=NA)
  lines(x=pred_surv1$u,y=pred_surv1$mean,col='blue',lwd=2)
  text(x=median(pred_surv1$u),y=0.99,labels =ifelse(t_0<=t_1,"Event 1","Event 2"), pos = 3, col = "black", cex = 1.6,font=2)
  axis(1, at = xticks)
  box(lwd=1.3)
  par(mar=c(0,0,0,4.5))
  xticks<-pretty(c(pred_surv2$u))
  plot(x=pred_surv2$u,y=pred_surv2$mean,col='black',type='l',las=1,yaxt='n',xaxt='n',xaxs='i',ylab='',ylim=c(0,1))
  polygon(c(pred_surv2$u,rev(pred_surv2$u)),c(pred_surv2$upper,rev(pred_surv2$lower)),col=rgb(0.6,0,0.6,0.4),border = NA)
  lines(x=pred_surv2$u,y=pred_surv2$mean,col='red',lwd=2)
  text(x=median(pred_surv2$u),y=0.99,labels =ifelse(t_0<=t_1,"Event 2","Event 1"), pos = 3, col = "black", cex = 1.6,font=2)
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
