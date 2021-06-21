#' Conducts simulation of abundance estimation for Huggins or HugFullHet models
#'
#' Generates simulation replicates of  double observer survey data with a known (marked) and unmarked portion of the population with
#' or without group sizes. Fits a sequence of models, and reports average abundance, average std error, confidence interval
#' coverage for each model and model average estimate.
#'
#' @param reps number of simulation replicates
#' @param n number of unmarked groups
#' @param nradios number of marked groups
#' @param beta vector of parameters for detection model specified in gen_formula
#' @param gen_formula formula specified as character for generating detection probabilities for simulated data
#' @param p0 proportion of population that has 0 probability of detection
#' @param data dataframe containing covariate data if something other than loc or cov (unmodelled heterogeneity); number of rows should match n+nradios
#' @param fit_formula character portion of the formula used to fit to observed double observer data excluding heterogeneity portion; use lngs for log group size
#' @param use_gs if TRUE, generate random group sizes from gamma distribution (integerized) with parameters in gs_param below
#' @param gs_param vector of 2 parameters for gamma distribution; gs_param[1]=scale and gs_param[2]=shape;see ?rgamma
#' @param which logical vector indicating which of the 6 models are fitted to the data
#' @param seed if NULL uses random seed, but if a value is used it can be used to repeat results from same seed
#' @return list containing: 1) abundance estimate (either groups or animals) , 2) std error, 3) lower 95 percent confidence limit,
#' 4) upper 95 percent confidence limit, 5) cover = 1 if TrueN in interval, 6) missed_low=1 if TrueN less than lower limit, and
#' 7) missed_high=1 if TrueN greater than upper limit.
#' @author Jeff Laake
#' @export
#' @examples
#' # runs 3 replicates as an example with parameters used to generate example_data and
#' # fits models with no terms other than loc and residual heterogeneity and others
#' # that added some combination of veg and lngs.
#'
#' simres=simhet(3,150,25,beta=c(-.5,-.5,-1,2,.5),data=data.frame(veg=runif(175)),
#' gen_formula="~-1+loc+veg+cov+lngs",fit_formula=c("veg+lngs","lngs","veg",""),
#' use_gs=TRUE,gs_param=c(1,10),seed=9234397)
simhet=function(reps,n,nradios,beta,gen_formula="~-1+loc+cov",p0=0,data=NULL,fit_formula="",use_gs=FALSE,gs_param=c(1,15),which=c(TRUE,TRUE,rep(FALSE,4)),seed=NULL)
{
  if(!is.null(seed)&is.numeric(seed))set.seed(seed)
  if(missing(reps))stop("must specify value for reps")
  if(missing(n))stop("must specify value for n - number of non-radio groups")
  if(missing(nradios))stop("must specify value for nradio - number of radio groups")
  if(missing(beta))stop("must specify value for beta - vector of parameters")
  if(use_gs&is.null(gs_param))stop("must specify value for gs_param if using group size (use_gs==TRUE")
  # if data is null use default of specifying location (loc) - 1/2; otherwise add to data
  if(is.null(data))
    data=data.frame(loc=factor(c(rep(1,n+nradios),rep(2,n+nradios))))
  else
  {
    savedata=data
    if(nrow(data)!=n+nradios)
      stop("\n number of rows in data must be n+nradios\n")
    else
    {
      data=rbind(data,data)
      data$loc=factor(c(rep(1,n+nradios),rep(2,n+nradios)))
    }
  }
  totaln=n+nradios
  TrueNgrp=totaln
  simresults=vector("list",length=reps)
  # loop over simulation reps
  for( i in 1:reps)
  {
    cat("\n i = ",i)
    # if use_gs, generate group sizes and add to data
    if(use_gs)
    {
      counts=ceiling(rgamma(nradios+n,shape=gs_param[1],scale=gs_param[2])) # generate group sizes
      counts_uncond=counts[1:nradios]
      counts_cond=counts[(nradios+1):totaln]
      TrueN=sum(counts)
      data$lngs=log(c(counts_uncond,counts_cond))
    } else
      TrueN=TrueNgrp
    # Include unmodelled covariate (-1,1) uniform variable in p but don't use in the model
    data$cov=rep(2*runif(totaln)-1,2)
    # Values of p for both observers
    pdesign=model.matrix(as.formula(gen_formula),data)
    if(length(beta)!=ncol(pdesign))
      stop(paste("length of beta vector ",length(beta), " does not match columns of design matrix ",ncol(pdesign)))
    p=plogis(pdesign%*%beta)
    #detection result for each observer
    ch1=rbinom(totaln,1,p[1:totaln]*(1-p0))
    ch2=rbinom(totaln,1,p[(totaln+1):(2*totaln)]*(1-p0))
    # unconditional (radio) and conditional (non-radio) capture histories
    ch_uncond=paste("1",ch1[1:nradios],ch2[1:nradios],sep="")
    ch_cond=paste("0",ch1[(nradios+1):totaln],ch2[(nradios+1):totaln],sep="")
    # throw away the "000" observations
    seen=ch_cond!="000"
    ch_cond=ch_cond[seen]
    # chi-square test of independence with unconditional observations -
    ch1=as.numeric(substr(ch_uncond[1:nradios],2,2))
    ch2=as.numeric(substr(ch_uncond[1:nradios],3,3))
    # chi sq test of independence
    xtab=table(ch1,ch2)
    suppressWarnings(Xsq<-chisq.test(xtab))
    # Use 3 characters with the first position always 1 to be equivalent to the
    # unconditional radio critters which is equivalent to a permanent mark in tag loss model;
    # having the 00 entries allows unmodelled heterogeneity to be estimated because
    # you get "100" capture histories. Without the "00" values you would have to assume
    # independence to estimate that product and unmodelled heterogeneity would cause bias.
    # Exactly equivalent to the tag loss problem. define an individual covariate het which
    # is 1 when first observer (position 2) detects the object (eg there is a 1 in position 2)
    # and 0 otherwise
    #
    # Create data frame
    if(!use_gs)
       df=data.frame(ch=c(ch_uncond,ch_cond),
                  het=sapply(strsplit(c(ch_uncond,ch_cond),""),function(x) as.numeric(x[2]=="1")),
                  stringsAsFactors = FALSE)
    else
    {
      df=data.frame(ch=c(ch_uncond,ch_cond), counts=c(counts_uncond,counts_cond[seen]),
                    het=sapply(strsplit(c(ch_uncond,ch_cond),""),function(x) as.numeric(x[2]=="1")),
                    stringsAsFactors = FALSE)
      df$lngs=log(df$counts)
    }
    df$type=factor(c(rep("marked",length(ch_uncond)),rep("unmarked",length(ch_cond))))
    df=cbind(df,rbind(savedata[1:nradios,,drop=FALSE],savedata[(nradios+1):totaln,,drop=FALSE][seen,,drop=FALSE]))

    # fit selection of models which returns marklist and abundance estimates
    results=fit_models(data=df,formulas=fit_formula,which=which,sim=TRUE)

    # compute coverage of confidence intervals
    coverage=function(N_lcl,N_ucl,TrueN)
    {
      return(data.frame(cover=as.numeric(TrueN>=N_lcl&TrueN<=N_ucl),
      missed_low=as.numeric(TrueN!=0)*as.numeric(TrueN<N_lcl),
      missed_high=as.numeric(TrueN>N_ucl)))

    }

    Nhat_group=cbind(results$Nhat_group,coverage(results$Nhat_group$lcl,results$Nhat_group$ucl,TrueNgrp))
    Nhat_group_modavg=cbind(results$Nhat_grp_modavg,coverage(results$Nhat_grp_modavg$lcl,results$Nhat_grp_modavg$ucl,TrueNgrp))
    if(use_gs)
    {
      Nhat=cbind(results$Nhat,coverage(results$Nhat$lcl,results$Nhat$ucl,TrueN))
      Nhat_modavg=cbind(results$Nhat_modavg,coverage(results$Nhat_modavg$lcl,results$Nhat_modavg$ucl,TrueN))
    }else
    {
      Nhat_modavg=NULL
      Nhat=NULL
    }
    #store index of best model
    best=as.numeric(row.names(results$marklist$model.table[1,]))

    # store simulation results for this rep
    simresults[[i]]=list(Xsq=Xsq,results=results$marklist,best=best,Nhatgrp=Nhat_group,Nhat=Nhat,TrueNgrp=TrueNgrp,TrueN=TrueN)
  }

  # summary of number of times each model was selected
  mod_sel=table(sapply(simresults,function(x)x$best))
  #percentage of cases chi-square rejected at 0.05 level
  chisq.p=mean(sapply(simresults, function(x)x$Xsq$p.value<0.05))*100
  #Average N of groups from simulations for each model
  models=simresults[[1]]$results$model.table
  modelnames=as.character(models[order(rownames(models)),"p"])
  Ngrp=data.frame(Model=1:nrow(simresults[[1]]$Nhatgrp),averageNgrp=rowMeans(sapply(simresults,function(x) x$Nhatgrp[order(x$Nhatgrp[,"Model#"]),"Nhat"])),
  se.averageNgrp=apply(sapply(simresults,function(x) x$Nhatgrp[order(x$Nhatgrp[,"Model#"]),"Nhat"]),1,function(x) sqrt(var(x)/reps)))
  Ngrp$PRB.Ngrp=100*(Ngrp$averageNgrp-TrueNgrp)/TrueNgrp
  Ngrp$coverage.Ngrp=rowMeans(sapply(simresults,function(x) x$Nhatgrp[order(x$Nhatgrp[,"Model#"]),"cover"]))
  Ngrp$missedlow.Ngrp=rowMeans(sapply(simresults,function(x) x$Nhatgrp[order(x$Nhatgrp[,"Model#"]),"missed_low"]))
  Ngrp$missedhigh.Ngrp=rowMeans(sapply(simresults,function(x) x$Nhatgrp[order(x$Nhatgrp[,"Model#"]),"missed_high"]))
  rownames(Ngrp)=c(modelnames,"Model average")
  if(use_gs)
  {
    #Average N from simulations for each model
    Nhat=data.frame(Model=1:nrow(simresults[[1]]$Nhat),averageNhat=rowMeans(sapply(simresults,function(x) x$Nhat[order(x$Nhat[,"Model#"]),"Nhat"])),
                    se.averageNhat=apply(sapply(simresults,function(x) x$Nhat[order(x$Nhat[,"Model#"]),"Nhat"]),1,function(x) sqrt(var(x)/reps)))
    Nhat$PRB.Nhat=100*(Nhat$averageNhat-TrueN)/TrueN
    Nhat$coverage.Nhat=rowMeans(sapply(simresults,function(x) x$Nhat[order(x$Nhat[,"Model#"]),"cover"]))
    Nhat$missedlow.Nhat=rowMeans(sapply(simresults,function(x) x$Nhat[order(x$Nhat[,"Model#"]),"missed_low"]))
    Nhat$missedhigh.Nhat=rowMeans(sapply(simresults,function(x) x$Nhat[order(x$Nhat[,"Model#"]),"missed_high"]))
    rownames(Nhat)=c(modelnames,"Model average")
   } else
  {
     Nhat=NULL
   }
  return(list(simresults=simresults,mod_sel=mod_sel,chisq.p=chisq.p,
              Ngrp=Ngrp,Nhat=Nhat))
}

