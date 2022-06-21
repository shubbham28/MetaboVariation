
#' @title
#' Lists flagged individuals in tabular format.
#' @description
#' A function to list the flagged individuals resulting from fitting a Bayesian generalised linear model (BGLM) model to repeated measures data.
#'
#' @param model An object of class \code{\link{MetaboVariation}}.
#' @param data A data frame containing data of all variables used in the BGLM along with a column detailing the individual ids.
#' @param individual_id A character string detailing the name of the column in the data frame that contains the individual ids.
#' @param threshold A positive integer Individuals that have been flagged in less than _threshold_ metabolites will not be included. By default, it is set to 2.
#'
#' @return The function will provide two types of lists based on the type of model object.
#' \itemize{
#'   \item Single metabolite model - If the model contains result for single metabolite, the function will return a list of flagged individuals along with their metabolite levels across different time points.
#'   \item Multiple metabolites model - If the model contain results for multiple metabolites, the function will return a list of individuals flagged in those metabolites more than _threshold_ times.
#' }
#' @seealso \code{\link{MetaboVariation}}
#' @export
#'
#' @examples
#' \dontrun{
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = getmetabolites(list = metabolite_list)
#' individual_id = "Individual_id"
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id
#' ,metabolite = metabolites,full_posterior = TRUE)
#' flaggedIndividuals(model,data = metabol.data,
#' individual_id = individual_id)
#' }
flaggedIndividuals <- function(model,data,individual_id,threshold = 2){
  if(!inherits(model,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
   #Single metabolite model
  if(inherits(model,"meta.single_model")){
    flagged = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=2),"identifier" = model$result[,6])
    flagged = as.factor(unique(flagged[flagged$identifier ==0,1]))
    flagged = data[data[[individual_id]] %in% flagged,c(individual_id,na.omit(stringr::str_extract(colnames(data),stringr::regex(paste(".*",model$metabolite,".*",sep = "")))))]
    rownames(flagged)=NULL
    return(as.data.frame(flagged))
  }

  # Multiple metabolite model
  else if(inherits(model,"meta.multi_model")) {
    individuals = c()
    for (met in names(model)) {
      status = cbind.data.frame(stringr::str_split_fixed(rownames(model[[met]]$result),pattern=" ",n=2),model[[met]]$result[,6])
      colnames(status) = c("id","timepoint","status")
      status = tidyr::spread(status,timepoint,status)
      individuals = cbind(individuals,apply(!status[,-1],2,function(x){ifelse(x==FALSE | is.na(x),0,1)}))
    }
    plotindividuals = c()
    timepoints = ncol(individuals)%/%length(names(model))
    timepoints = 1:timepoints
    for (met in names(model)) {
      plotindividuals = cbind(plotindividuals,apply(individuals[,paste(model[[met]]$metabolite,"_",timepoints,sep = "")],1,sum))
    }
    rownames(plotindividuals) = status$id
    summ = apply(individuals, 1,sum)
    plotindividuals = plotindividuals[summ >= threshold,]
    output = cbind(rownames(plotindividuals),data.frame(plotindividuals, row.names=NULL))
    colnames(output) = c("Individual_id",names(model))
    rownames(output)=NULL
    return(as.data.frame(output))
  }


}


#Throw


