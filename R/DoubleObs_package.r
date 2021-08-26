#' Introduction
#'
#' DoubleObs provides two primary capabilities: 1) analyze hybrid double observer survey data in which some
#' of the animal groups are marked with a radio/satellite transmitter so it is known whether they were seen or not,
#' and 2) to simulate data from a hybrid double-observer survey and fit models to explore aspects of this type of analysis.
#' We describe these two capabilities below.
#'
#' Analysis: Model Fitting/Abundance Estimation
#' Data
#'   There are a few requirements for the dataframe used in the analysis. Each record in the dataframe is for a single
#'   group of animals. The following 3 fields are required with specific names for each record:
#'    1) ch - must be a character string of length 3 with values of 0 or 1 in each position. The first character is 1 for groups that
#'    are marked (eg. an individual in group has a transmitter) and 0 otherwise.  The second character is a 1 if the front observer(s)
#'    detected the group and 0 otherwise.  The third character is the same for the rear observer(s).
#'    2) count - the count of animals in the group
#'    2) type - a factor variable with values "marked" and "unmarked".
#'   Any number of additional variables can be defined for use in the model fitting (eg lngs = log(count), vegetation type/cover etc).
#' Model Fitting
#'   Models are fitted with MARK using the RMark interface with the Huggins model primarily. As such, various steps need to be taken to use RMark and an
#'   understanding of RMark is helpful but not entirely necessary because the DoubleObs package contains functions: prep_data, prep_ddl and
#'   fit_models as wrappers for the RMark code. See the help files ?prep_data, ?prep_ddl and ?fit_models for more details. The prep_data
#'   and prep_ddl run the functions process.data and make.design.data in RMark and add some specific fields that are useful in fitting models to
#'   hybrid double observer survey data. Alternatively you can use the RMark functions separately and you can replace fit_models with calls to mark.wrapper and
#'   mark to fit models directly. The example ?hybrid_analysis shows different ways of fitting models.
#' Abundance estimation
#'   Two functions Nhat_group_marklist and Nhat_marklist in this package will construct abundance of groups and animals, respectively,
#'   for a set of models (marklist) and will compute the model averaged estimates as well.
#' While this package is mostly intended for double observer hybrid models, we have provided examples for a double observer survey without
#' marked groups ?double_observer_analysis and a sightability analysis in which a sample of marked animals is used to create the detection
#' model and those predictions are used for the survey data ?sightability_analysis.
#'
#' Simulation:
#' The function simhet provides a simulation capability to generate hybrid double observer survey data with residual
#' heterogeneity. The function senerates simulation replicates of  double observer survey data with a known (marked) and unmarked portion of the population with
#' or without group sizes. Fits a sequence of models, and reports average abundance, average std error, confidence interval
#' coverage for each model and model average estimate.
#'
#' @name DoubleObs
#' @author Jeff Laake
NULL

#' Example Double Observer Hybrid Data Analysis
#'
#' An example of double observer data with 25 marked (radio) groups and 150 other groups (no radio)
#' in the population. Data were generated with the following code using simhet function.
#'
#'simres=simhet(1,150,25,beta=c(-.5,-.5,-1,2,.5),data=data.frame(veg=runif(175)),
#'              gen_formula="~-1+loc+veg+cov+lngs",use_gs=TRUE,gs_param=c(1,10),seed=9234397)
#'example_data=simres[[1]][[1]]$results[[1]]$data$data
#'
#'
#'The groups sizes were randomly generated from a gamma distribution with scale=10 and shape=1 (gs_param).
#'The true simulated population size is 175 groups and 1732 animals. A uniform random veg variable, and
#'unmodelled heterogeneity variable cov (uniform -1,1) and natural log of group size are used to
#'generate p for each location(loc).
#'
#'@name hybrid_analysis
#'@aliases example_data
#'@docType data
#'@format A data frame with 129 observations of 7 variables:
#' \describe{\item{ch}{a character vector containing the three character encounter history; first occasion indicates whether the group was marked or not}
#' \item{count}{number of animals in the group}
#' \item{het}{variable used to model heterogeneity (re-capture probability.}
#' \item{lngs}{ natural log of count}
#' \item{veg}{random vegetation variable}
#' \item{type}{"unmarked" or "marked"}
#' }
#' @keywords datasets
#' @examples
#' data(example_data)
#' # example using fit_models
#' fit_models(data=example_data,formula=c("veg+lngs","lngs","veg",""))
#' #example using RMark directly
#' dp=prep_data(example_data,model="Huggins")
#' # variables het, rear and loc are created in prep_ddl which call make.design.data in RMark
#' # and then creates other variables used to fit the models
#' ddl=prep_ddl(dp)
#' fit_do=function()
#' {
#'    # tag loss approach
#'    p.1=list(formula=~veg+lngs+loc+c:het:rear,share=TRUE)
#'    p.2=list(formula=~lngs+loc+c:het:rear,share=TRUE)
#'    p.3=list(formula=~veg+loc+c:het:rear,share=TRUE)
#'    p.4=list(formula=~loc+c:het:rear,share=TRUE)
#'    # hybrid approach
#'    p.5=list(formula=~veg+lngs+loc+type,share=TRUE)
#'    p.6=list(formula=~lngs+loc+type,share=TRUE)
#'    p.7=list(formula=~veg+loc+type,share=TRUE)
#'    p.8=list(formula=~loc+type,share=TRUE)
#'    cml=RMark::create.model.list("Huggins")
#'    results=RMark::mark.wrapper(cml,data=dp,ddl=ddl)
#'    return(results)
#' }
#' results=fit_do()
#' results
#' # Group abundance
#' Nhat_group_marklist(results,example_data)
#' # Group abundance estimate from MARK for best model
#' colSums(results[[1]]$results$derived[[1]])
#' # Group abundance not using 100 set (marked groups missed)
#' # This can be useful if not all marked missed groups can be found after survey
#' Nhat_group_marklist(results,example_data,marked_known=FALSE)
#' # MARK doesn't consider group size so their is no equivalent for total abundance
#' Nhat_marklist(results,example_data)
#' # Animal abundance not using 100 set (marked groups missed)
#' # This can be useful if not all marked missed groups can be found after survey
#' Nhat_marklist(results,example_data,marked_known=FALSE)
NULL


#' Example Double Observer Data Analysis
#'
#' An example of double observer data analysis but not using marked individuals so residual heterogeneity
#' cannot be considered.
#' Data were generated with the following code using simhet function.
#'
#'simres=simhet(1,150,25,beta=c(-.5,-.5,-1,2,.5),data=data.frame(veg=runif(175)),
#'              gen_formula="~-1+loc+veg+cov+lngs",use_gs=TRUE,gs_param=c(1,10),seed=9234397)
#'example_data=simres[[1]][[1]]$results[[1]]$data$data
#'
#'
#'The groups sizes were randomly generated from a gamma distribution with scale=10 and shape=1 (gs_param).
#'The true simulated population size is 175 groups and 1732 animals. A uniform random veg variable, and
#'unmodelled heterogeneity variable cov (uniform -1,1) and natural log of group size are used to
#'generate p for each location(loc).
#'
#'@name double_observer_analysis
#'@docType data
#'@format A data frame with 129 observations of 7 variables:
#' \describe{\item{ch}{a character vector containing the three character encounter history; first occasion indicaes whther the group was marked or not}
#' \item{count}{number of animals in the group}
#' \item{het}{variable used to model heterogeneity (re-capture probability.}
#' \item{lngs}{ natural log of count}
#' \item{veg}{random vegetation variable}
#' \item{type}{"unmarked" or "marked"}
#' }
#' @keywords datasets
#' @examples
#' data(example_data)
#' df=example_data
#' # Only use second and third capture history positions - first position was for radio(marked)
#' df$ch=substr(df$ch,2,3)
#' # Ignore 00 values: these were marked groups that were not detected
#' df=df[df$ch!="00",]
#' # all are considered unmarked but need to create empty marked group with call to factor
#' df$type="unmarked"
#' df$type=factor(df$type,levels=c("marked","unmarked"))
#' # Use process.data in RMark but use groups="type" and allgroups=TRUE to create group structure
#' # even though it is not used. In addition, begin.time=2. These options make double observer
#' # analysis work with abundance estimation code Nhat functions.
#' dp=RMark::process.data(df,model="Huggins",groups="type",allgroups=TRUE,begin.time=2)
#' # Create default design data
#' ddl=RMark::make.design.data(dp)
#' # Define create rear 0/1 variable - time==3 (second occasion)
#' ddl$p$rear=ifelse(ddl$p$time==2,0,1)
#' ddl$c$rear=1
#' # show design data
#' ddl$p
#' ddl$c
#' # Create function to fit a sequence of models; must assume independence in this model
#' # so cannot include c:het:rear or type in models to account for residual heterogeneity
#' fit.models=function()
#' {
#'    p.1=list(formula=~veg+lngs,share=TRUE)
#'    p.2=list(formula=~lngs+rear,share=TRUE)
#'    p.3=list(formula=~veg+rear,share=TRUE)
#'    p.4=list(formula=~veg+lngs+rear,share=TRUE)
#'    p.5=list(formula=~rear,share=TRUE)
#'    cml=RMark::create.model.list("Huggins")
#'    results=RMark::mark.wrapper(cml,data=dp,ddl=ddl)
#'    return(results)
#' }
#' results=fit.models()
#' results
#' # Group abundance
#' Nhat_group_marklist(results,df)
#' # Group abundance estimate from MARK for best model
#' results[[1]]$results$derived[[1]][2,]
#' # MARK doesn't consider group size so there is no equivalent for total abundance
#' Nhat_marklist(results,df)
NULL



#' Example Sightability Data Analysis
#'
#' An example of sightability data analysis using marked individuals to build a logistic regression model
#' and then estimating abundance from those seen.
#'
#' Data were generated with the following code using simhet function.
#'
#'simres=simhet(1,150,25,beta=c(-.5,-.5,-1,2,.5),data=data.frame(veg=runif(175)),
#'              gen_formula="~-1+loc+veg+cov+lngs",use_gs=TRUE,gs_param=c(1,10),seed=9234397)
#'example_data=simres[[1]][[1]]$results[[1]]$data$data
#'
#'
#'The groups sizes were randomly generated from a gamma distribution with scale=10 and shape=1 (gs_param).
#'The true simulated population size is 175 groups and 1732 animals. A uniform random veg variable, and
#'unmodelled heterogeneity variable cov (uniform -1,1) and natural log of group size are used to
#'generate p for each location(loc).
#'
#'@name sightability_analysis
#'@docType data
#'@format A data frame with 129 observations of 7 variables:
#' \describe{\item{ch}{a character vector containing the three character encounter history; first occasion indicaes whther the group was marked or not}
#' \item{count}{number of animals in the group}
#' \item{het}{variable used to model heterogeneity (re-capture probability.}
#' \item{lngs}{ natural log of count}
#' \item{veg}{random vegetation variable}
#' \item{type}{"unmarked" or "marked"}
#' }
#' @keywords datasets
#' @examples
#' data(example_data)
#' df=example_data
#' # Only use "sum" of second and third capture history positions; first position for radio(marked)
#' df$seen=ifelse((substr(df$ch,2,2)=="1") | (substr(df$ch,3,3)=="1"),1,0)
#' # For logistic model development only use marked groups.
#' lrdata=df[df$type=="marked",]
#' # restrict df to only those seen for prediction
#' df=df[df$seen==1,]
#' # Only fit a single model for this example
#' lrmodel=glm(seen~veg+lngs,family="binomial",data=lrdata)
#' # get predicted value of p and std error for each survey observation
#' pred=predict(lrmodel,df,type="response",se.fit=TRUE)
#' pdot=pred$fit
#' # construct variance-covariance matrix of predictions p
#' vcov=summary.glm(lrmodel)$cov.unscaled
#' dm=model.matrix(~veg+lngs,data=df)
#' deriv=dm*pdot*(1-pdot)
#' vc_pdot=deriv%*%vcov%*%t(deriv)
#' # construct variance
#' deriv=matrix(-1/pdot^2,nrow=1)
#' #' # Group abundance
#' Ngrp=sum(1/pdot)
#' se_Ngrp=sqrt(sum((1-pdot)/pdot^2)+deriv%*%vc_pdot%*%t(deriv))
#' Mt1=nrow(df)
#' Ngrp=Ngrp-Mt1
#' cv_Ngrp=se_Ngrp/Ngrp
#' C=exp(1.96*sqrt(log(1+cv_Ngrp^2)))
#' Ngrp_lcl=Ngrp/C+Mt1
#' Ngrp_ucl=Ngrp*C+Mt1
#' Ngrp=Ngrp+Mt1
#' paste("Ngrp = ",round(Ngrp), " se = ", round(se_Ngrp), " CI = (",round(Ngrp_lcl),",
#' ",round(Ngrp_ucl),")",sep="")
#'
#' # Animal abundance
#' deriv=matrix(-df$count/pdot^2,nrow=1)
#' N=sum(df$count/pdot)
#' se_N=sqrt(sum(df$count^2*(1-pdot)/pdot^2)+deriv%*%vc_pdot%*%t(deriv))
#' Mt1=sum(df$count)
#' N=N-Mt1
#' cv_N=se_N/N
#' C=exp(1.96*sqrt(log(1+cv_N^2)))
#' N_lcl=N/C+Mt1
#' N_ucl=N*C+Mt1
#' N=N+Mt1
#' paste("N = ",round(N), " se = ", round(se_N), " CI = (",round(N_lcl),",",round(N_ucl),")",sep="")
NULL



