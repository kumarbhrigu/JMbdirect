#' Function for bootstrapped confidence interval
#'
#' @param object fitted model
#' @param future_time time sequence at which estimates are required
#'
#' @return Returns bootstraped confidence interval for model fitted through FastJM
#' @import jmBIG
#' @export
#'
bootciJMCS<-function(object,future_time){
  if(!inherits(object,"survfitJMCS"))
    stop("\n Not a 'survfitJMCS' object.\n")
  longdata<-object$others$jmcs_others$others$dtlong
  survdata<-object$others$jmcs_others$others$dtsurv
  id<-object$others$jmcs_others$others$id
  idNumber<-object$others$ids
  jmcsModel<-object$others$jmcs_others
  obs.time<-object$others$obs.time
  longm<-object$others$jmcs_others$others$longm
  survm<-object$others$jmcs_others$others$survm
  rd<-object$others$jmcs_others$others$rd
  samplesize<-object$others$jmcs_others$samplesize
  lastid<-max(survdata[[id]])
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
  bootstrapped_data <- bootstrap_longitudinal_survival(longdata,
                                                       survdata,
                                                       n_bootstrap = 10,id=id,
                                                       idNumber = idNumber)

  fit2<-list();P2<-list();P3<-list()
  for(i in 1:10){
    P2<-try(survfitJMCS(model=jmcsModel,ids=idNumber,
                        u=future_time,
                        obs.time=obs.time),silent=TRUE)
    fit2[[i]]<-try(jmcsBig(dtlong=data.frame(bootstrapped_data[[i]]$longitudinal),
                           dtsurv = data.frame(bootstrapped_data[[i]]$survival),
                           longm=longm,survm=survm,
                           rd= rd,samplesize=200,id=id),silent=TRUE)

    P2[[i]]<-try(survfitJMCS(model=fit2[[i]],ids=c(lastid),
                             u=future_time,
                             obs.time=obs.time),silent=TRUE)
    P3[[i]]<-try(P2[[i]]$P1$Pred[[as.character(lastid)]][,2],silent=TRUE)

  }

  CIdata<-Reduce('cbind',P3)
  qCIdata<-data.frame(Times=P2[[i]]$P1$Pred[[as.character(lastid)]]$times,LL=apply(CIdata,1,function(x){quantile(x,0.025)}),Med=apply(CIdata,1,function(x){quantile(x,0.5)}),UL=apply(CIdata,1,function(x){quantile(x,0.975)}))
  result<-list()
  result$lastid<-lastid
  result$bootCI<-qCIdata
  result$P2<-P2
  result$P3<-P3
  class(result)<-"cisurvfitJMCS"
  result

}
