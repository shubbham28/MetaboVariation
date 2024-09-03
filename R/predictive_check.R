globalVariables(c("Model","Timepoint","MAD"))
#' @title Use posterior predictive checks to assess the model fit of a MetaboVariation model.
#' @description
#' This function performs posterior predictive checks to assess the fit of a MetaboVariation model by generating replicate datasets from the posterior predictive distribution and computing the mean absolute differences between correlation matrices for the original and replicate datasets at each time point.
#' @param data A data frame containing all variables to be used in the BGLM.
#' @param model An object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' @param timepoints A vector of integer detailing the time points present in the data. For example. if there are four timepoints present in the data, the user should specify c(1,2,3,4).
#' @param metabolites A string or vector of strings containing the names of metabolites to model.
#' @param replication Number of replicate datasets to generate. Defaults to 100. Should be less than effective sample size, calculated as (number of iterations - number of warmup iterations) / thinning rate as declared in class \code{\link{MetaboVariation}}.
#'
#' @return The function returns two values.
#' \itemize{
#' \item plot - A box plot displaying the mean absolute differences between correlation matrices of the original data and the replicate data sets.
#' \item mad - A data frame containing mad values for each time point for each replicate data set.
#' }
#'  along with the MAD values for each replicate data set and each timepoint.
#' @export
#'
#' @examples
#'\dontrun{
#' # Load the simulated data and extract the metabolites names.
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#'
#' # Run MetaboVariation on first three metabolites.
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3],type="dependent")
#'
#' # Perform posterior predictive check
#' check = predictive_check(metabol.data, model, c(1,2,3,4), metabolites)
#' # To plot the box plot
#' check$plot
#' # To see MAD values for the model.
#' check$mad
#' }
predictive_check <- function(data,model,timepoints,metabolites,replication=100){
  if(!inherits(model,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }

  original = list()
  data = stats::na.omit(data)
  for (i in timepoints){
    original[[i]] = stats::cor(data[,colnames(data) %like% paste0("_",i)])
  }

  mad_model1 = matrix(0,nrow = replication,ncol = length(timepoints))
  replicated_model1 = prediction(model$model[[1]])$posterior[sample(nrow(model$model[[1]]$Sol),replication,replace = FALSE),]
  for (i in 2:length(model$model)) {
    replicated_model1 = replicated_model1 + prediction(model$model[[i]])$posterior[sample(nrow(model$model[[i]]$Sol),replication,replace = FALSE),]
  }
  replicated_model1 = replicated_model1 / length(model$model)

  names = rownames(model$result)
  names = names[order(stringr::str_sub(names,-3))]
  colnames(replicated_model1) = names


  for (r in 1:replication) {

  sub_model1 = replicated_model1[r,]
  sub_model1 = matrix(sub_model1,ncol=length(metabolites))
  colnames(sub_model1) = metabolites
  rownames(sub_model1) = unique(stringr::str_sub(names,end=-5))

  for (i in timepoints){
    mad_model1[r,i] = mean(abs(stats::cor(sub_model1[rownames(sub_model1) %like% paste0(" ",i),]) - original[[i]]))
  }
  }

  mad_model1 = as.data.frame(mad_model1)
  colnames(mad_model1) = timepoints
  mad_model1 = reshape2::melt(mad_model1,variable.name = "Timepoint",value.name = "MAD")

 plotdf = mad_model1
  plot = ggplot2::ggplot(plotdf,ggplot2::aes(x=Timepoint,y=MAD,color = Timepoint))+ggplot2::geom_boxplot() + ggplot2::labs(x="Timepoint",y="MAD") + ggplot2::theme_classic()+
    ggplot2::theme(legend.position = "none")
  return(list(plot = plot, mad = plotdf))
}
