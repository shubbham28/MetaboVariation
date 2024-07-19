`%>%` <- magrittr::`%>%`
`%like%` <- data.table::`%like%`
globalVariables(c("metabolite","flag","original","timepoint","lower","upper","nonzero_count","n","individual","total_count","percentage","value"))

#' @title Heatmap of individuals flagged across metabolites.
#' @description
#' Visualise a heatmap that shows individuals flagged in metabolite across time points.
#'
#' @param model An object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' @param interval The interval for the highest posterior distribution (HPD) that you want to visualise.
#' @param threshold A positive integer indicating the minimum number of metabolites an individual must be flagged across to be included in the heatmap.
#'
#' @return Returns a heatmap that illustrates the number of individuals flagged in all the metabolites across time points.
#' @export
#'
#' @examples

#' \dontrun{
#' # Load the simulated data and extract the metabolites names.
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#'
#' # Run the MetaboVariation.
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], covariates = covariates,cutoff=c(0.95,0.975,0.99))
#'
#' # Visualise the heatmap.
#' metabolite.heatmap(model,interval = 0.975,threshold=50)
#'}


metabolite.heatmap <- function(model,interval = NULL, threshold = 1){
  if(!inherits(model,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
  if(model$type!="dependent"){
    stop("Function is only used with dependent model")
  }
  test= model$result
  subdf = test[,paste0("flag",interval)]
  subdf = data.frame(subdf)
  colnames(subdf) = c("flag")
  subdf$individual = gsub("(.*) .* .*","\\1",rownames(subdf))
  subdf$timepoint = gsub(".* (.*) .*","\\1",rownames(subdf))
  subdf$metabolite = gsub(".* .* (.*)","\\1",rownames(subdf))
  #subdf$flag = sapply(subdf$flag,function(x){ifelse(x==0,1,0)})
  timepoints=length(unique(subdf$timepoint))
  metabolite_flagged <- stats::aggregate(flag ~ individual + metabolite, data = subdf, FUN = sum)
  #metabolite_flagged$flag = sapply(metabolite_flagged$flag,function(x){ifelse(x>=(timepoints*as.numeric(type)),"Yes","No")})
  individuals = unique(metabolite_flagged[metabolite_flagged$flag>=(timepoints*0.01),"individual"])
  metabolite_flagged$flag = factor(metabolite_flagged$flag)
  metabolite_flagged = metabolite_flagged[metabolite_flagged$individual %in% individuals,]
  # flagged_metabolites <- metabolite_flagged %>%
  #   dplyr::group_by(metabolite) %>%
  #   dplyr::summarize(nonzero_count = sum(flag != 0), total_count = dplyr::n()) %>%
  #   dplyr::mutate(percentage = nonzero_count / total_count) %>%
  #   dplyr::filter(percentage > 0.1) %>%
  #   dplyr::distinct(metabolite)
#  metabolite_flagged = metabolite_flagged[metabolite_flagged$metabolite %in% flagged_metabolites,]
  flagged_individuals <- metabolite_flagged %>%
    dplyr::group_by(individual) %>%
    dplyr::summarize(nonzero_count = sum(flag != 0), total_count = dplyr::n()) %>%
    # dplyr:: mutate(percentage = nonzero_count / total_count) %>%
    dplyr::filter(nonzero_count >= threshold) %>%
    dplyr::distinct(individual)
   metabolite_flagged = metabolite_flagged[metabolite_flagged$individual %in% unlist(flagged_individuals),]


  color <- grDevices::colorRampPalette(c("white", "#2B7CE9"))(length(levels(metabolite_flagged$flag))+1)
  names(color) = c("",levels(metabolite_flagged$flag))
  ggplot2::ggplot(metabolite_flagged,ggplot2::aes(x=metabolite,y=individual,fill=flag))+ggplot2::geom_tile(color = "white")+
  ggplot2::scale_fill_manual(values = color)+ggplot2::theme_minimal()+ggplot2::theme(text = ggplot2::element_text(size = 20),axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),legend.position = "top")  + ggplot2::labs(x = "Metabolites",y="Flagged individuals",fill = "Number of flagged time points")
}
