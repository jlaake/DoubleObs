#' Prepares data and design data (ddl)
#'
#' prep_data: Certain variable names and values are hard coded and this prepares the data and ddl to make sure
#' they are setup correctly and fixes specific real values in the ddl.
#'
#' 1) The data needs to contain a variable het which is constructed if not available in data,
#' 2) A factor variable type with values "marked" and "unmarked" must be in data to define groups in data,
#' 3) The encounter history named ch must be of type character with 3 characters each,
#' 4) The model argument must be either Huggins or HugFullHet, and
#' 5) The extra ... arguments are passed to process.data but cannot contains "groups" which is pre-specified.
#'
#' prep_ddl: Prepares the ddl to make sure and fixes specific real values in the ddl.
#' 1) for p and c parameters creates variable loc (for location) with values "Radio","Front" and "Rear" for 3 characters in ch,
#' 2) for p and c parameters, creates a numeric variable "rear" which is 1 for loc=="Rear" and 0 otherwise,
#' 3) creates variable fix for parameters p and c, to fix real parameter p to 1 for radio for type=="marked" and 0 for type="unmarked",
#' 4) for HugFullHet model (pi in ddl), the pi for type=="unmarked" is set to 0 because that mixture has p=0 so by definition can't be in unmarked which only contain those seen,
#' 5) for HugFullHet model (pi in ddl), set p=0 (except for Radio loc) for mixture because they cannot be seen by definition
#'
#' @usage
#'   prep_data(data,model,...)
#'   prep_ddl(dp,...)
#'
#' @param data dataframe containing double observer data with field ch (capture history), freq (optional), type ="u" (marked) or "c" (unmarked), count if more than one animal in each sighting. This is passed to process.data after checks.
#' @param model either "Huggins" or "HugFullHet"
#' @param dp processed data passed to make.design.data after checks
#' @param ... extra parameters passed to either process.data or make.design.data in RMark
#' @export prep_data
#' @export prep_ddl
#' @aliases prep_ddl prep_data
#' @return processed data list from prep_data and design data list (ddl) from prep_ddl
#'
prep_data=function(data,model,...)
{
  # check fields in data
  if(is.null(data$ch)) stop("\nmissing ch (capture history field)\n")
  if(!is.character(data$ch)) stop("\nch (capture history field) must be type character\n")
  if(any(nchar(data$ch)!=3))stop("each ch string must be of length 3 characters")
  if(is.null(data$type)) stop("\nmissing type field to indicate marked and unmarked groups\n")
  if(!is.factor(data$type))
  {
    message("type not a factor variable. changing to a factor variable")
    data$type=factor(data$type)
  }
  if(any(!c("marked","unmarked")%in%levels(data$type)) | length(levels(data$type))!=2 )
    stop("\ntype variable must have only 2 levels marked and unmarked")
  # create het variable
  data$het=sapply(strsplit(data$ch,""),function(x) as.numeric(x[2]=="1"))
  # check model; must be either Huggins or HugFullHet
  if(!model%in%c("Huggins","HugFullHet"))stop("\nmodel must be either Huggins or HugFullHet")
  # check to make sure groups argument is not in ...
  if("groups" %in% names(list(...))) stop("\ndo not specify groups argument\n")
  # Process data with HugFullHet model with 2 mixtures with probabilities pi and 1-pi
  if(model=="HugFullHet")
     dp <- process.data(data, model = "HugFullHet", mixtures=2, groups = c("type"),...)
  else
     dp <- process.data(data, model = "Huggins", groups = c("type"),...)
  return(dp)
}
prep_ddl=function(dp,...)
{
  # create default design data
  ddl <- make.design.data(dp,...)
  #create loc from Time
  loc=factor(c(ddl$p$Time,ddl$c$Time+1),labels=c("Radio","Front","Rear"))
  # change the initial level for the loc field to "front" because first value is fixed
  loc=relevel(loc,"Front")
  ddl$p$loc=loc[1:nrow(ddl$p)]
  ddl$c$loc=loc[(nrow(ddl$p)+1):length(loc)]
  # create a 0/1 variable (rear) that is 1 for rear location (time==3) and 0 otherwise; this is done
  # for p and c because we will share the parameters in the same columns of the design
  # matrix for these 2 parameters. To do that the design dataframes have to be row binded
  # together and must have the same columns to do so
  ddl$p$rear=ifelse(ddl$p$time==3,1,0)
  ddl$c$rear=ifelse(ddl$c$time==3,1,0)
  # create a field call fix which is NA when the real parameter should be estimated;
  # otherwise it is the value you want to specify. Here we are using it to fix p=1 for
  # first position (radio location; time==1) for radio groups and p=0 for non-radio groups
  ddl$p$fix=NA
  ddl$p$fix[ddl$p$time==1&ddl$p$type=="marked"]=1
  ddl$p$fix[ddl$p$time==1&ddl$p$type=="unmarked"]=0
  ddl$c$fix=NA
  #If HugFullHet model
  if("pi"%in%names(ddl))
  {
    # fix pi=0 for unmarked groups as none can be seen by default (p=0)
    ddl$pi$fix=NA
    ddl$pi$fix[ddl$pi$type=="unmarked"]=0
    # fix p=0 for first mixture because they are the proportion of the population with p=0
    ddl$p$fix[ddl$p$mixture==1&ddl$p$loc!="Radio"]=0
    ddl$c$fix[ddl$c$mixture==1&ddl$c$loc!="Radio"]=0
  }
  return(ddl)
}
