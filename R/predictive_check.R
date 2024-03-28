globalVariables(c("Model","variable"))
#' Performs posterior predictive check for two MetaboVariation objects
#'
#' @param data A data frame containing data of all variables to be used in the BGLM.
#' @param model1 An object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' @param model2 A different object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' @param timepoints A list of time points present in the data.
#' @param metabolites A string or vector of string containing the name of metabolites to model.
#' @param replication number of replicated data set to make. Defaults to 100. Should be less than effective sample size.
#'
#' @return A boxplot comparing the differences between two \code{\link{MetaboVariation}} objects. Replicated datasets are created from both models,
#' generating a correlation matrix of the original dataset and replicated dataset for each timepoint. The boxplot displays the mean absolute difference for each timepoint.
#' @export
#'
#' @examples
#'\dontrun{
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1], covariates = covariates)
#' model1 = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3],type="dependent")
#' model2 = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3],type="independent")
#' predictive_check(metabol.data, model1, model2, c(1,2,3,4), metabolites)
#' }
predictive_check <- function(data,model1,model2,timepoints,metabolites,replication=100){
  if(!inherits(model1,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
  if(!inherits(model2,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
  original = list()
  data = stats::na.omit(data)
  for (i in timepoints){
    original[[i]] = stats::cor(data[,colnames(data) %like% paste0("_",i)])
  }

  mad_model1 = matrix(0,nrow = replication,ncol = length(timepoints))
  mad_model2 = matrix(0,nrow = replication,ncol = length(timepoints))


  replicated_model1 = prediction(model1$model[[1]])$posterior[sample(nrow(model1$model[[1]]$Sol),replication,replace = FALSE),]
  for (i in 2:length(model1$model)) {
    replicated_model1 = replicated_model1 + prediction(model1$model[[i]])$posterior[sample(nrow(model1$model[[i]]$Sol),replication,replace = FALSE),]
  }
  replicated_model1 = replicated_model1 / length(model1$model)

  replicated_model2 = prediction(model2$model[[1]])$posterior[sample(nrow(model2$model[[1]]$Sol),replication,replace = FALSE),]
  for (i in 2:length(model2$model)) {
    replicated_model2 = replicated_model2 + prediction(model2$model[[i]])$posterior[sample(nrow(model2$model[[i]]$Sol),replication,replace = FALSE),]
  }
  replicated_model2 = replicated_model2 / length(model2$model)
  names = rownames(model1$result)
  names = names[order(stringr::str_sub(names,-3))]
  colnames(replicated_model1) = names
  colnames(replicated_model2) = names

  for (r in 1:replication) {

  sub_model1 = replicated_model1[r,]
  sub_model1 = matrix(sub_model1,ncol=length(metabolites))
  colnames(sub_model1) = metabolites
  rownames(sub_model1) = unique(stringr::str_sub(names,end=-5))

  for (i in timepoints){
    mad_model1[r,i] = mean(abs(stats::cor(sub_model1[rownames(sub_model1) %like% paste0(" ",i),]) - original[[i]]))
  }

  sub_model2 = replicated_model2[r,]
  sub_model2 = matrix(sub_model2,ncol=length(metabolites))
  colnames(sub_model2) = metabolites
  rownames(sub_model2) = unique(stringr::str_sub(names,end=-5))

  for (i in timepoints){
    mad_model2[r,i] = mean(abs(stats::cor(sub_model2[rownames(sub_model2) %like% paste0(" ",i),]) - original[[i]]))
  }
  }

  mad_model1 = as.data.frame(mad_model1)
  mad_model2 = as.data.frame(mad_model2)
  colnames(mad_model1) = timepoints
  colnames(mad_model2) = timepoints
  mad_model1 = reshape2::melt(mad_model1)
  mad_model2 = reshape2::melt((mad_model2))
  plotdf = rbind(mad_model1,mad_model2)
  plotdf$Model = c(rep("Model 1",nrow(mad_model1)),rep("Model 2",nrow(mad_model2)))
  ggplot2::ggplot(plotdf,ggplot2::aes(x=variable,y=value,color = Model))+ggplot2::geom_boxplot() + ggplot2::labs(colour ="Model",x="Timepoint",y="MAD") + ggplot2::theme_classic()

}
