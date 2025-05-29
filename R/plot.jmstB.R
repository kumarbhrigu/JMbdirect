#' @title Prediction plot from \code{jmstB()}
#' @param x fitted model object
#' @param y newdata
#' @param ... others
#' @note
#' In the example code we use newdata as the data for ID 2 in the PBC2 dataset, it has follow up information till
#' 8.832. Now suppose we want to look at the survival of ID 2 under joint model
#' 1 after time 4 and for joint model 2 after time 9. For that we created the
#' newdata as if the individual is followed till for a time period
#' less than min(4,9).
#' @return Returns prediction plot for the newdata using the model fitted through \code{jmstB()}.
#' @import rstanarm
#' @import jmBIG
#' @importFrom stats median
#' @examples
#'  \donttest{
#' ##
#' library(JMbayes2)
#' library(rstanarm)
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
#' model_jmstBdirect<-jmstB(
#'   dtlong=pbc2_new,
#'   dtsurv = new_pbc2id,
#'   longm=list(serBilir ~ drug * year+(year|id),albumin~drug+year+(year|id)),
#'   survm=list(Surv(years,status2) ~ drug,Surv(time_2,status_2) ~ drug),
#'   timeVar="year",
#'   id='id',
#'   refresh=400,
#'   nchain=1)
#' t0<-4
#' nd<-pbc2[pbc2$id %in% c(2), ]
#' nd<-nd[nd$year<t0,]
#' nd$status2<-0
#' nd$years<-t0
#' nd$time_2<-9
#' nd$status_2<-0
#' plot(x=model_jmstBdirect,y = nd)
#' ##
#' }
#' @rdname plot.jmstB
#' @method plot jmstB
#' @export
plot.jmstB<-function(x,y,...){
  if(!inherits(x,'jmstB'))
    stop("\n Not a 'jmbBdirect' object.\n")
  result<-x
  model1<-x$model1
  model2<-x$model2
  y<-as.data.frame(y)
  newdata<-y
  IDvar<-x$IDvar
  ID<-as.numeric(unique(y[,IDvar]))
  timeVar<-x$timeVar

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
  oldpar <- par(no.readonly = TRUE)
  K<-2
  m<-cbind(1:K,rep(K+1,K),rep(K+2,K))
  xticks<-pretty(c(newdata1$year,newdata2$year))
  widths<-c(0.3,0.35,0.35)
  layout(m,widths=widths)
  par(mar=c(0,4.5,0,0),oma=c(4,0,3,0),cex.axis=1.1,font.axis=2,font.lab=2,font.main=2)
  plot(x=newdata1$year,y=newdata1$yvalue1,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=newdata1$year,y=newdata1$yvalue1,pch=8,col='red')
  mtext(first_long,line=2.5,side=2,font=2)
  box(lwd=1.3)
  par(mar = c(0, 4.5, 0, 0))
  plot(x=newdata2$year,y=newdata2$yvalue2,type='l',las=1,ylab='',xaxt='n',lwd=1.5)
  points(x=newdata2$year,y=newdata2$yvalue2,pch=8,col='red')
  #points(x=newdata2$year,y=newdata2$yvalue2,pch=5,col='green')
  mtext(scnd_long,line=2.5,side=2,font=2)
  axis(1, at = xticks)
  box(lwd=1.3)
  par(mar=c(0,0,0,0))
  xticks<-pretty(c(psurv1$year))
  plot(x=psurv1$year,y=psurv1$survpred,col='black',type='l',las=1,yaxt='n',xaxt='n',xaxs='i',ylab='',ylim=c(0,1))
  polygon(c(psurv1$year,rev(psurv1$year)),c(psurv1$ci_ub,rev(psurv1$ci_lb)),col=rgb(0,0,0.6,0.4),border=NA)
  lines(x=psurv1$year,y=psurv1$survpred,col='blue',lwd=2)
  text(x=median(psurv1$year),y=0.99,labels = ifelse(t_0<=t_1,"Event 1","Event 2"), pos = 3, col = "black", cex = 1.6,font=2)
  axis(1, at = xticks)
  box(lwd=1.3)
  par(mar=c(0,0,0,4.5))
  xticks<-pretty(c(psurv2$year))
  plot(x=psurv2$year,y=psurv2$survpred,col='black',type='l',las=1,yaxt='n',xaxt='n',xaxs='i',ylab='',ylim=c(0,1))
  polygon(c(psurv2$year,rev(psurv2$year)),c(psurv2$ci_ub,rev(psurv2$ci_lb)),col=rgb(0.6,0,0.6,0.4),border = NA)
  lines(x=psurv2$year,y=psurv2$survpred,col='red',lwd=2)
  text(x=median(psurv2$year),y=0.99,labels = ifelse(t_0<=t_1,"Event 2","Event 1"), pos = 3, col = "black", cex = 1.6,font=2)
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


utils::globalVariables(c('yfit','yvalue1','ci_lb','ci_ub','survpred'))



