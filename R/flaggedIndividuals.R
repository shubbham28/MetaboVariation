
#' @title
#' Gives flagged individuals in a tablular format
#' @description
#' From the class "MetaboVariation", the results contain the individuals that have been flagged as they have values outside the 95% central credible interval based on the results.
#' This function provides the details of those individuals in a tabular format.
#'
#' @param model An object of class "MetaboVariation"
#' @param data An object of class data.frame (or one that can be coerced to that class) containing data of all variables used for modelling.
#' @param individual_id A character string that denotes the column name of individual ids present in the dataset passed.
#' @param covariates A list of column names present in the data that denotes the covariates data for the individuals in the data.
#' @param threshold A number to filter out individuals that have been flagged in less than "threshold" metabolites. By default, it is set at 2.
#'
#' @return A data.frame of the flagged individuals that have values outside of 95% credible interval.
#' @export
#'
#' @examples
#' \dontrun{
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Sample.Id"
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id
#' ,metabolite = metabolites[1:3],covariates = covariates,full_posterior = TRUE)
#' flaggedIndividuals(model,data = metabol.data,
#' individual_id = individual_id,covariates=covariates)
#' }
flaggedIndividuals <- function(model,data,individual_id,covariates=NULL,threshold = 2){
  if(!inherits(model,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
   #Single metabolite model
  if(inherits(model,"meta.single_model")){
    flagged = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=2),"identifier" = model$result[,6])
    flagged = as.factor(unique(flagged[flagged$identifier ==0,1]))
    flagged = data[data[[individual_id]] %in% flagged,c(individual_id,covariates,na.omit(stringr::str_extract(colnames(data),stringr::regex(paste(".*",model$metabolite,".*",sep = "")))))]
    return(as.data.frame(flagged))
  }

  # Multiple metabolite model
  else if(inherits(model,"meta.multi_model")) {

    individuals = c()
    for (met in names(model)) {
      status = cbind.data.frame(stringr::str_split_fixed(rownames(model[[met]]$result),pattern=" ",n=2),model[[met]]$result[,6])
      colnames(status) = c("id","timepoint","status")
      time_data = as.numeric(gsub("[^0-9]", "",unique(status$timepoint)))
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
    individuals = data[data[[individual_id]] %in% as.factor(rownames(plotindividuals)),c(individual_id,covariates)]
    return(as.data.frame(individuals))
  }


}
