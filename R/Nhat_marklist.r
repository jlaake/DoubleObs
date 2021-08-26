#' Abundance estimation for Huggins or HugFullHet models from marklist of fitted models
#'
#' Uses (Nhat_group) and/or animal (Nhat) abundance to estimate abundance for a list
#' of fitted mark models. Works with results of fit_models
#'
#'
#' @usage
#' Nhat_group_marklist(marklist,data,marked_known=TRUE)
#' Nhat_marklist(marklist,data,marked_known=TRUE)
#' mod_avgN(marklist,Nmatrix,Ngrp=TRUE,data)
#'
#' @aliases Nhat_group_marklist Nhat_marklist mod_avgN
#' @param marklist list of fitted mark models
#' @param data dataframe containing data used to fit model or subset of that data
#' @param marked_known if FALSE includes marked groups into estimation,  otherwise adds number in marked groups to abundance
#' @param Ngrp if TRUE, model average is for group abundance
#' @param Nmatrix matrix of estimates in order of model.table (return value from Nhat_group_marklist and Mhat_marklist)
#'
#' @return matrix of group abundance estimates (Nhat_group_marklist) or animal abundance (Nhat_marklist) with std errors, and confidence intervals for each model and model averaged values
#' @author Jeff Laake
#' @export Nhat_marklist Nhat_group_marklist mod_avgN
#'
Nhat_marklist=function(marklist,data,marked_known=TRUE)
{
  if(any(sapply(marklist[1:(length(marklist)-1)],function(x) ncol(x$data$group.covariates)!=1))) stop("\nOnly works with one group variable named 'type'\n")
  if(any(sapply(marklist[1:(length(marklist)-1)],function(x) is.null(x$data$data$type)))) stop("\n group variable named 'type' not found in models\n")
  if(any(sapply(marklist[1:(length(marklist)-1)],function(x) !x$model%in%c("Huggins","HugFullHet"))))stop("\nOnly works with Huggins or HugFullHet models\n")
  if(is.null(data$type)) stop("\nmissing type field to indicate marked and unmarked groups\n")
  if(!is.factor(data$type))
  {
    message("type not a factor variable. changing to a factor variable")
    data$type=factor(data$type)
  }
  if(any(!c("marked","unmarked")%in%levels(data$type)) | length(levels(data$type))!=2 )
    stop("\ntype variable must have only 2 levels marked and unmarked")
  if(!marked_known)
  {
    data=data[!(data$type=="marked"&substr(data$ch,2,3)=="00"),]
    data$type="unmarked"
  }
  known=sum(data$count[data$type=="marked"])
  if(class(marklist)!="marklist") stop("\nmarklist must contain a list of mark models and model.table as created by mark.wrapper\n")
  Nmatrix=matrix(NA,nrow=nrow(marklist$model.table),ncol=4)
  # Np0 logical vector in order of models in marklist as to whether adjustment for pi should be done (TRUE)
  #Nhybrid  logical vector in order of models in marklist as to whether the hybrid approach is used for abundance (TRUE - uses type-=="u" parameters)
  Np0=rep(FALSE,nrow(marklist$model.table))
  Nhybrid=rep(FALSE,nrow(marklist$model.table))
  for(i in 1:(length(marklist)-1))
  {
    if(is.null(marklist[[i]]$parameters$pi$fixed)&marklist[[i]]$model!="Huggins") Np0[i]=TRUE
    if("type" %in% all.vars(marklist[[i]]$parameters$p$formula)) Nhybrid[i]=TRUE
  }
  for(i in 1:length(Np0))
  {
    ddl=marklist[[i]]$design.data
    if(Np0[i])
    {
      indices=c(ddl$pi$model.index[ddl$pi$type=="marked"],ddl$p$model.index[ddl$p$type=="unmarked"&ddl$p$mixture==2&ddl$p$time%in%2:3])
      if(Nhybrid[i])stop("\nModels with estimated pi should not contain type in formula for p\n")
      Nest=Nhat0(marklist[[i]],indices=indices,data=data,known=known)
    }else
    {
      if(Nhybrid[i])
      {
        if(marklist[[i]]$model=="HugFullHet")
           indices=c(ddl$p$model.index[ddl$p$type=="marked"&ddl$p$mixture==2&ddl$p$time%in%2:3])
        else
          indices=c(ddl$p$model.index[ddl$p$type=="marked"&ddl$p$time%in%2:3])
        Nest=Nhat(marklist[[i]],indices=indices,data=data,known=known)
      }
      else
      {
        if(marklist[[i]]$model=="HugFullHet")
           indices=c(ddl$p$model.index[ddl$p$type=="unmarked"&ddl$p$mixture==2&ddl$p$time%in%2:3])
        else
          indices=c(ddl$p$model.index[ddl$p$type=="unmarked"&ddl$p$time%in%2:3])
        Nest=Nhat(marklist[[i]],indices=indices,data=data,known=known)
      }
    }
    Nmatrix[i,]=c(Nest[[1]][1],as.numeric(Nest[[2]]),as.numeric(Nest[[3]]),as.numeric(Nest[[4]]))
  }
  model_order=as.numeric(rownames(marklist$model.table))
  Nmatrix=Nmatrix[model_order,]
  Nmatrix=cbind(model_order,Nmatrix)
  colnames(Nmatrix)=c("Model#","Nhat","se","lcl","ucl")
  rownames(Nmatrix)=paste(c("","p0 ")[as.numeric(Np0)+1],marklist$model.table$p,sep="")
  Nhat_modavg=mod_avgN(marklist,Nmatrix,FALSE,data=data)
  Nmatrix=rbind(Nmatrix,cbind("Model#"=nrow(Nmatrix)+1, Nhat_modavg))
  rownames(Nmatrix)[nrow(Nmatrix)]="Model average"
  return(Nmatrix)
}

Nhat_group_marklist=function(marklist,data,marked_known=TRUE)
{
  if(any(sapply(marklist[1:(length(marklist)-1)],function(x) ncol(x$data$group.covariates)!=1))) stop("\nOnly works with one group variable named 'type'\n")
  if(any(sapply(marklist[1:(length(marklist)-1)],function(x) is.null(x$data$data$type)))) stop("\n group variable named 'type' not found in models\n")
  if(any(sapply(marklist[1:(length(marklist)-1)],function(x) !x$model%in%c("Huggins","HugFullHet"))))stop("\nOnly works with Huggins or HugFullHet models\n")
  if(is.null(data$type)) stop("\nmissing type field to indicate marked and unmarked groups\n")
  if(!is.factor(data$type))
  {
    message("type not a factor variable. changing to a factor variable")
    data$type=factor(data$type)
  }
  if(any(!c("marked","unmarked")%in%levels(data$type)) | length(levels(data$type))!=2 )
    stop("\ntype variable must have only 2 levels marked and unmarked")
  if(!marked_known)
  {
    data=data[!(data$type=="marked"&substr(data$ch,2,3)=="00"),]
    data$type="unmarked"
  }
  if(class(marklist)!="marklist") stop("\nmarklist must contain a list of mark models and model.table as created by mark.wrapper\n")
  known=nrow(data[data$type=="marked",])
  Nmatrix=matrix(NA,nrow=nrow(marklist$model.table),ncol=4)
  # Np0 logical vector in order of models in marklist as to whether adjustment for pi should be done (TRUE)
  #Nhybrid  logical vector in order of models in marklist as to whether the hybrid approach is used for abundance (TRUE - uses type-=="u" parameters)
  Np0=rep(FALSE,nrow(marklist$model.table))
  Nhybrid=rep(FALSE,nrow(marklist$model.table))
  for(i in 1:(length(marklist)-1))
  {
    if(is.null(marklist[[i]]$parameters$pi$fixed)&marklist[[i]]$model!="Huggins") Np0[i]=TRUE
    if("type" %in% all.vars(marklist[[i]]$parameters$p$formula)) Nhybrid[i]=TRUE
  }
  for(i in 1:length(Np0))
  {
    ddl=marklist[[i]]$design.data
    if(Np0[i])
    {
      indices=c(ddl$pi$model.index[ddl$pi$type=="marked"],ddl$p$model.index[ddl$p$type=="unmarked"&ddl$p$mixture==2&ddl$p$time%in%2:3])
      if(Nhybrid[i])stop("\nModels with estimated pi should not contain type in formula for p\n")
      Nest=Nhat_grp0(marklist[[i]],indices=indices,data=data,known=known)
    }else
    {
      if(Nhybrid[i])
      {
        if(marklist[[i]]$model=="HugFullHet")
          indices=c(ddl$p$model.index[ddl$p$type=="marked"&ddl$p$mixture==2&ddl$p$time%in%2:3])
        else
          indices=c(ddl$p$model.index[ddl$p$type=="marked"&ddl$p$time%in%2:3])
        Nest=Nhat_grp(marklist[[i]],indices=indices,data=data,known=known)
      }
      else
      {
        if(marklist[[i]]$model=="HugFullHet")
          indices=c(ddl$p$model.index[ddl$p$type=="unmarked"&ddl$p$mixture==2&ddl$p$time%in%2:3])
        else
          indices=c(ddl$p$model.index[ddl$p$type=="unmarked"&ddl$p$time%in%2:3])
        Nest=Nhat_grp(marklist[[i]],indices=indices,data=data,known=known)
      }
    }
    Nmatrix[i,]=c(Nest[[1]][1],as.numeric(Nest[[2]]),as.numeric(Nest[[3]]),as.numeric(Nest[[4]]))
  }
  model_order=as.numeric(rownames(marklist$model.table))
  Nmatrix=Nmatrix[model_order,]
  Nmatrix=cbind(model_order,Nmatrix)
  colnames(Nmatrix)=c("Model#","Nhat","se","lcl","ucl")
  rownames(Nmatrix)=paste(c("","p0 ")[as.numeric(Np0)+1],marklist$model.table$p,sep="")
  Nhat_modavg=mod_avgN(marklist,Nmatrix,data=data)
  Nmatrix=rbind(Nmatrix,cbind("Model#"=nrow(Nmatrix)+1, Nhat_modavg))
  rownames(Nmatrix)[nrow(Nmatrix)]="Model average"
  return(Nmatrix)
}

mod_avgN=function(marklist,Nmatrix,Ngrp=TRUE,data)
{
  if(Ngrp)
    Mt1=nrow(data[data$type=="unmarked",])
  else
    Mt1=sum(data$count[data$type=="unmarked"])
  if(Ngrp)
    known=nrow(data[data$type=="marked",])
  else
    known=sum(data$count[data$type=="marked"])
  # drop any models with negative variances
  drop=rep(FALSE,nrow(marklist$model.table))
  rownums=as.numeric(rownames(marklist$model.table))
  for(i in 1:length(rownums))
    drop[i]=any(diag(marklist[[rownums[i]]]$results$beta.vcv)<0)
  #adjust weights after dropping any
  marklist$model.table$weight[drop]=0
  marklist$model.table$weight=marklist$model.table$weight/sum(marklist$model.table$weight)
  Nmatrix[drop,"se"]=0
  # compute model average after dropping any models
  modelavgN=sum((Nmatrix[,"Nhat"]-known)*marklist$model.table$weight)
  deviations=Nmatrix[,"Nhat"]-known-modelavgN
  se.modelavgN=sum(marklist$model.table$weight*sqrt(deviations^2+Nmatrix[,"se"]^2))
  cvsq.modavgN=(se.modelavgN/(modelavgN-Mt1))^2
  C=exp(1.96*sqrt(log(1+cvsq.modavgN)))
  N_lcl=(modelavgN-Mt1)/C+Mt1+known
  N_ucl=(modelavgN-Mt1)*C+Mt1+known
  return(data.frame(Nhat=modelavgN+known,se=se.modelavgN,lcl=N_lcl,ucl=N_ucl))
}
