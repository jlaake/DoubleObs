#' Abundance estimation for Huggins or HugFullHet models
#'
#' Computes estimate of group (Nhat_group) or animal (Nhat) abundance for Huggins model or
#' for HugFullHet with 2 mixtures in which pi is the proportion of animals with p=0.
#'
#'
#' @usage
#' Nhat_grp(mod,indices,data,known)
#'
#' Nhat(mod,indices,data,known)
#'
#' Nhat0(mod,indices,data,known)
#'
#' Nhat_grp0(mod,indices,data,known)
#'
#' @aliases Nhat_grp Nhat Nhat0 Nhat_grp0
#' @param mod fitted Huggins or HugFullHet model for double observer data
#' @param indices model.index values of real parameters for p and pi(for Nhat0 and Nhat_grp0) from design data
#' @param data dataframe containing data used to fit model
#' @param known known marked animals (e.g. radio collar); either number of marked groups or sum of animals in marked groups
#' @return list containing: 1) abundance estimate (either groups or animals) , 2) std error, 3) lower 95 percent confidence limit,
#' 4) upper 95 percent confidence limit.
#' @author Jeff Laake
#' @importFrom RMark covariate.predictions
# function to compute group abundance, std error and 95% log-normal confidence interval for Huggins estimator
Nhat_grp=function(mod,indices,data,known)
{
  # get real estimates of p1 and p2 for each non-radio observation (type=="unmarked")
  # p for observer 1 and 2 are produced for each conditional observation; the indices vary
  # because the hybrid observer model uses the p1 and p2 from the radio observations with the
  # non-radio observations
  vars=all.vars(mod$model.parameters$p$formula)
  vars=c("type",vars[vars%in%names(data)])
  xx=covariate.predictions(mod,data=data[data$type=="unmarked",vars,drop=FALSE],indices=indices)
  # split estimates into pairs and then use sapply to compute pdot for each observation
  pdot=sapply(split(xx$estimates$estimate,factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(x[1]+x[2]-x[1]*x[2]))
  #abundance of groups
  Ngrp=sum(1/pdot)

  # the variance code uses the chain-rule by first computing the var-cov matrix of the pdot estimates (vc_pdot) and then from it computing the
  # the variance of sum(1/pdot) for the n values (70 in this case)
  nr=nrow(xx$vcv)
  # first step
  deriv=sapply(split(xx$estimates$estimate,factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(c(1-x[2],1-x[1],rep(0,nr))))
  deriv=matrix(as.vector(deriv),byrow=T,ncol=nr)[-(nr/2+1),]
  vc_pdot=deriv%*%xx$vcv%*%t(deriv)
  #second step - first term is variance given estimates and second term due to estimation ofN'N parameters
  deriv=matrix(-1/pdot^2,nrow=1)
  se_Ngrp=sqrt(sum((1-pdot)/pdot^2)+deriv%*%vc_pdot%*%t(deriv))

  # confidence interval with f0 approach
  Mt1=nrow(data[data$type=="unmarked",])
  Ngrp=Ngrp-Mt1
  cv_Ngrp=se_Ngrp/Ngrp
  C=exp(1.96*sqrt(log(1+cv_Ngrp^2)))
  Ngrp_lcl=Ngrp/C+Mt1+known
  Ngrp_ucl=Ngrp*C+Mt1+known
  Ngrp=Ngrp+Mt1+known
  return(list(Ngrp,se_Ngrp,Ngrp_lcl=Ngrp_lcl,Ngrp_ucl=Ngrp_ucl))
}

# function to compute animal abundance, std error and 95% log-normal confidence interval for Huggins estimator
Nhat=function(mod,indices,data,known)
{
  if(is.null(data$count))stop("\ndata must contain a field called count to estimate animal abundance\n")
  # get real estimates of p1 and p2 for each non-radio observation (type=="unmarked")
  # p for observer 1 and 2 are produced for each conditional observation; the indices vary
  # because the hybrid observer model uses the p1 and p2 from the radio observations with the
  # non-radio observations
  vars=all.vars(mod$model.parameters$p$formula)
  vars=c("type",vars[vars%in%names(data)])
  xx=covariate.predictions(mod,data=data[data$type=="unmarked",vars,drop=FALSE],indices=indices)
  # split estimates into pairs and then use sapply to compute pdot for each observation
  pdot=sapply(split(xx$estimates$estimate,factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(x[1]+x[2]-x[1]*x[2]))
  #abundance of animals
  N=sum(data$count[data$type=="unmarked"]/pdot)

  # the variance code uses the chain-rule by first computing the var-cov matrix of the pdot estimates (vc_pdot) and then from it computing the
  # the variance of sum(1/pdot) for the n values (70 in this case)
  nr=nrow(xx$vcv)
  # first step
  deriv=sapply(split(xx$estimates$estimate,factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(c(1-x[2],1-x[1],rep(0,nr))))
  deriv=matrix(as.vector(deriv),byrow=T,ncol=nr)[-(nr/2+1),]
  vc_pdot=deriv%*%xx$vcv%*%t(deriv)
  #second step - first term is variance given estimates and second term due to estimation ofN'N parameters
  deriv=matrix(-data$count[data$type=="unmarked"]/pdot^2,nrow=1)
  se_N=sqrt(sum(data$count[data$type=="unmarked"]^2*(1-pdot)/pdot^2)+deriv%*%vc_pdot%*%t(deriv))

  # confidence interval with f0 approach
  Mt1=sum(data$count[data$type=="unmarked"])
  N=N-Mt1
  cv_N=se_N/N
  C=exp(1.96*sqrt(log(1+cv_N^2)))
  N_lcl=N/C+Mt1+known
  N_ucl=N*C+Mt1+known
  N=N+Mt1+known
  return(list(N=N,se_N=se_N,N_lcl=N_lcl,N_ucl=N_ucl))
}

# function to compute group abundance, std error and 95% log-normal confidence interval for HugFullHet with
# 2 mixtures and mixture 1 has p=0
Nhat_grp0=function(mod,indices,data,known)
{
  # get real estimates of p1 and p2 for each non-radio observation (type=="unmarked")
  # p for observer 1 and 2 are produced for each conditional observation; the indices vary
  # because the hybrid observer model uses the p1 and p2 from the radio observations with the
  # non-radio observations
  vars=all.vars(mod$model.parameters$p$formula)
  vars=c("type",vars[vars%in%names(data)])
  df=data[data$type=="unmarked",vars,drop=FALSE]
  df=df[rep(1:nrow(df),each=2),,drop=FALSE]
  df$index=indices[2:3]
  df=rbind(df[1,,drop=FALSE],df)
  df$index[1]=indices[1]
  xx=covariate.predictions(mod,data=df)
  # split estimates into pairs and then use sapply to compute pdot for each observation
  pi=xx$estimates$estimate[1]
  pdot=(1-pi)*sapply(split(xx$estimates$estimate[-1],factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(x[1]+x[2]-x[1]*x[2]))
  #abundance of groups
  Ngrp=sum(1/pdot)

  # the variance code uses the chain-rule by first computing the var-cov matrix of the pdot estimates (vc_pdot) and then from it computing the
  # the variance of sum(1/pdot) for the n values (70 in this case)
  nr=nrow(xx$vcv)-1
  # first step
  deriv=sapply(split(xx$estimates$estimate[-1],factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(c(1-x[2],1-x[1],rep(0,nr))))
  deriv=matrix(as.vector(deriv),byrow=T,ncol=nr)[-(nr/2+1),]*(1-pi)
  deriv=cbind(-pdot,deriv)
  vc_pdot=deriv%*%xx$vcv%*%t(deriv)
  #second step - first term is variance given estimates and second term due to estimation ofN'N parameters
  deriv=matrix(-1/pdot^2,nrow=1)
  se_Ngrp=sqrt(sum((1-pdot)/pdot^2)+deriv%*%vc_pdot%*%t(deriv))

  # confidence interval with f0 approach
  Mt1=nrow(data[data$type=="unmarked",])
  Ngrp=Ngrp-Mt1
  cv_Ngrp=se_Ngrp/Ngrp
  C=exp(1.96*sqrt(log(1+cv_Ngrp^2)))
  Ngrp_lcl=Ngrp/C+Mt1+known
  Ngrp_ucl=Ngrp*C+Mt1+known
  Ngrp=Ngrp+Mt1+known
  return(list(Ngrp,se_Ngrp,Ngrp_lcl=Ngrp_lcl,Ngrp_ucl=Ngrp_ucl))
}

# function to compute animal abundance, std error and 95% log-normal confidence interval for HugFullHet with
# 2 mixtures and mixture 1 has p=0
Nhat0=function(mod,indices,data,known)
{
  if(is.null(data$count))stop("\ndata must contain a field called count to estimate animal abundance\n")
  # get real estimates of pi and p1 and p2 for each non-radio observation (type=="unmarked")
  # p for observer 1 and 2 are produced for each conditional observation; the indices vary
  # because the hybrid observer model uses the p1 and p2 from the radio observations with the
  # non-radio observations
  vars=all.vars(mod$model.parameters$p$formula)
  vars=c("type",vars[vars%in%names(data)])
  df=data[data$type=="unmarked",vars,drop=FALSE]
  df=df[rep(1:nrow(df),each=2),,drop=FALSE]
  df$index=indices[2:3]
  df=rbind(df[1,,drop=FALSE],df)
  df$index[1]=indices[1]
  xx=covariate.predictions(mod,data=df)
  # split estimates into pairs and then use sapply to compute pdot for each observation
  pi=xx$estimates$estimate[1]
  pdot=(1-pi)*sapply(split(xx$estimates$estimate[-1],factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(x[1]+x[2]-x[1]*x[2]))
  #abundance of groups
  N=sum(data$count[data$type=="unmarked"]/pdot)

  # the variance code uses the chain-rule by first computing the var-cov matrix of the pdot estimates (vc_pdot) and then from it computing the
  # the variance of sum(1/pdot) for the n values (70 in this case)
  nr=nrow(xx$vcv)-1
  # first step
  deriv=sapply(split(xx$estimates$estimate[-1],factor(rep(1:nrow(data[data$type=="unmarked",]),each=2))),function(x) return(c(1-x[2],1-x[1],rep(0,nr))))
  deriv=matrix(as.vector(deriv),byrow=T,ncol=nr)[-(nr/2+1),]*(1-pi)
  deriv=cbind(-pdot,deriv)
  vc_pdot=deriv%*%xx$vcv%*%t(deriv)
  #second step - first term is variance given estimates and second term due to estimation ofN'N parameters
  deriv=matrix(-data$count[data$type=="unmarked"]/pdot^2,nrow=1)
  se_N=sqrt(sum(data$count[data$type=="unmarked"]^2*(1-pdot)/pdot^2)+deriv%*%vc_pdot%*%t(deriv))

  # confidence interval with f0 approach
  Mt1=nrow(data[data$type=="unmarked",])
  N=N-Mt1
  cv_N=se_N/N
  C=exp(1.96*sqrt(log(1+cv_N^2)))
  N_lcl=N/C+Mt1+known
  N_ucl=N*C+Mt1+known
  N=N+Mt1+known
  return(list(N,se_N,N_lcl=N_lcl,N_ucl=N_ucl))
}

formatN_output=function(Nhat)
  paste("N = ",round(Nhat[[1]][1]), " se = ", round(as.numeric(Nhat[[2]])), " CI = (",round(as.numeric(Nhat[[3]])),",",round(as.numeric(Nhat[[4]])),")",sep="")
