#' Example Double Observer Data
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
#'unmodelled heterogeneity variable cov (uniform -1,1) and natural log of grop size are used to
#'generate p for each location(loc).
#'
#'@name example_data
#'@docType data
#'@format A data frame with 129 observations of 7 variables:
#' \describe{\item{ch}{a character vector containing the three character encounter history; first occasion indicaes whther the group was marked or not}
#' \item{counts}{number of animals in the group}
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
#' Nhat_group_marklist(results,example_data)
#' Nhat_marklist(results,example_data)
NULL

