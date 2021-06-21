#' Fits a sequence of heterogeneity HugFullHet models
#'
#' Fits a sequence of residual heterogeneity models for each formula specified. The models in order are the
#' following: 1) ~-1+loc+c:het:rear - tag loss model, 2) ~loc+type - original hybrid observer model,
#' 3) ~loc+type+c:het:rear - combined tagloss-hybrid approach, 4) ~loc - no residual heterogeneity (fitIndep==TRUE),
#' 5) ~loc+type+c:het:rear - tag loss with estimated pi - proportion with p=0, and 6) ~loc independence except with estimated pi - proportion with p=0
#' Each model is fitted via HugFullHet but if pi=0 this becomes the Huggins model. Only models 5 and 6 have an estimated pi in which
#' p=0 and for 1-pi mixture the p's are estimated. You can specify which of these 6 models are used with the logical argument "which".
#'
#' If you were to use formulas=c("","lngs") the function will fit 2*sum(which) models.
#'
#' Some of the heterogeneity models use the variable het in the data.  If it is not in the data, it is created from ch.
#'
#' @param data dataframe containing double observer data with field ch (capture history), freq (optional), type ="u" (marked) or "c" (unmarked), count if more than one animal in each sighting.
#' @param formulas character vector of formulas other than residual heterogeneity portion of the formula
#' @param which logical vector indicating which of the 6 models are fitted to the data
#' @param sim if TRUE, deletes output files
#' @return list containing: 1) marklist of fitted models, 2) Nhat_group, group abundance for each model, 3) Nhat, animal abundance if data$count exists.
#' @author Jeff Laake
#' @export
#' @importFrom RMark process.data make.design.data create.model.list mark.wrapper merge.mark
#' @importFrom stats as.formula chisq.test model.matrix plogis rbinom relevel rgamma runif var
fit_models=function(data,formulas="",which=c(TRUE,TRUE,rep(FALSE,4)),sim=FALSE)
{
  # prepare data
  dp=prep_data(data,model="HugFullHet")
  ddl=prep_ddl(dp)
  fit_non0=function(dp,ddl,formulas=NULL,which)
  {
    if(sum(which[1:4])==0)return(NULL)
    j=1
    for(i in 1:length(formulas))
    {
      if(formulas[i]!="")
      {
        if(which[1])
        {
          assign(paste("p",j,sep="."),list(formula=as.formula(paste("~-1+loc+c:het:rear",formulas[i],sep="+")),share=TRUE))
          j=j+1
        }
        if(which[2])
        {
          assign(paste("p",j,sep="."),list(formula=as.formula(paste("~loc+type",formulas[i],sep="+")),share=TRUE))
          j=j+1
        }
        if(which[3])
        {
          assign(paste("p",j,sep="."),list(formula=as.formula(paste("~loc+type+c:het:rear",formulas[i],sep="+")),share=TRUE))
          j=j+1
        }
        if(which[4])
        {
          assign(paste("p",j,sep="."),list(formula=as.formula(paste("~loc",formulas[i],sep="+")),share=TRUE))
          j=j+1
        }
      }else
      {
        if(which[1])
        {
          assign(paste("p",j,sep="."),list(formula=~-1+loc+c:het:rear,share=TRUE))
          j=j+1
        }
        if(which[2])
        {
          assign(paste("p",j,sep="."),list(formula=~-1+loc+type,share=TRUE))
          j=j+1
        }
        if(which[3])
        {
          assign(paste("p",j,sep="."),list(formula=~loc+type+c:het:rear,share=TRUE))
          j=j+1
        }
        if(which[4])
        {
          assign(paste("p",j,sep="."),list(formula=~loc,share=TRUE))
          j=j+1
        }
      }
    }
    pi.1=list(formula=~1,fixed=0)
    cml=create.model.list("HugFullHet")
    # setting delete=TRUE removes the external files that are not used
    if(!sim)
      results=mark.wrapper(cml,data=dp,ddl=ddl,brief=TRUE,delete=FALSE)
    else
      results=mark.wrapper(cml,data=dp,ddl=ddl,brief=TRUE,delete=TRUE)
    return(results)
  }
  fit0=function(dp,ddl,formulas=NULL,which)
  {
    if(sum(which[5:6])==0)return(NULL)
    j=1
    for(i in 1:length(formulas))
    {
      if(formulas[i]!="")
      {
        if(which[5])
        {
          assign(paste("p",j,"0",sep="."),list(formula=as.formula(paste("~-1+loc+c:het:rear",formulas[i],sep="+")),share=TRUE))
          j=j+1
        }
        if(which[6])
        {
          assign(paste("p",j,"0",sep="."),list(formula=as.formula(paste("~loc",formulas[i],sep="+")),share=TRUE))
          j=j+1
        }
      } else
      {
        if(which[5])
        {
          assign(paste("p",j,"0",sep="."),list(formula=~-1+loc+c:het:rear,share=TRUE))  #mod
          j=j+1
        }
        if(which[6])
        {
          assign(paste("p",j,"0",sep="."),list(formula=~loc,share=TRUE)) # independence model
          j=j+1
        }
      }
    }
    pi.1=list(formula=~1)
    cml=create.model.list("HugFullHet")
    # setting delete=TRUE removes the external files that are not used
    if(!sim)
      results=mark.wrapper(cml,data=dp,ddl=ddl,brief=TRUE,delete=FALSE)
    else
    results=mark.wrapper(cml,data=dp,ddl=ddl,brief=TRUE,delete=TRUE)
    return(results)
  }
  # get results for models without mixture
  results_non0=fit_non0(dp=dp,ddl=ddl,formulas,which=which)
  # get results for models with mixture
  results0=fit0(dp=dp,ddl=ddl,formulas,which=which)
  # merge model results
  if(is.null(results0))
    if(is.null(results_non0))
      stop("\nno models chosen\n")
    else
      results=results_non0
  else
    if(is.null(results_non0))
       results=results0
    else
       results=merge.mark(results_non0,results0)
  # Compute group abundance (Nhat_group) and animal abundance (Nhat),if data$count is found
  Nhat_group=Nhat_group_marklist(results,data=data)
  if(!is.null(data$count))
    Nhat=Nhat_marklist(results,data=data)
   else
    Nhat=NULL
  return(list(marklist=results,Nhat_group=Nhat_group,Nhat=Nhat))
}


