globalVariables(c("flagged_metabolites",".data"))

#' @title
#' Pair-wise metabolite heatmap.
#' @description
#' Plots a heatmap that shows the number of flagged individuals in which pairs of metabolites contribute to flagging the individuals across all time points.
#'
#' @param model An object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' @param interval.width The width of the highest posterior density (HPD) interval considered. Must be a numeric between 0 and 1 with default value of 0.95.
#'
#' @return Returns a heatmap that illustrates the number of individuals flagged in the pairs of metabolites across all time points.
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
#' # Run MetaboVariation on first three metabolites.
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], covariates = covariates,interval.width=c(0.95,0.975,0.99))
#'
#' # Visualise the metabolite pair heatmap.
#' metpair.heatmap(model,interval.width=0.975)
#'}

metpair.heatmap <-function(model,interval.width = 0.95){

  model = as.data.frame(model$result[])
  model$timepoint = gsub(".* (.*) .*","\\1",rownames(model))
  model$individual = gsub("(.*) .* .*","\\1",rownames(model))
  model$metabolite = gsub(".* .* (.*)","\\1",rownames(model))
  metabolites = unique(model$metabolite)
  flagged_dep = model %>%  dplyr::filter(.data[[paste0("flag",interval.width)]] == 1) %>%   dplyr::group_by(individual,timepoint) %>%  dplyr::summarise(flagged_metabolites = dplyr::n_distinct(metabolite),
                                                                                                                                                  metabolites_list = toString(sort(unique(metabolite))))

  flagged_dep$metabolites <- strsplit(as.character(flagged_dep$metabolites_list), ", ")
  binary_matrix <- sapply(metabolites, function(met) {
    sapply(flagged_dep$metabolites, function(x) as.integer(met %in% x))
  })

  co_occurrence_matrix <- t(binary_matrix) %*% binary_matrix
  co_occurrence_df <- as.data.frame(as.table(as.matrix(co_occurrence_matrix)))
  dist_matrix <- stats::dist(co_occurrence_matrix)
  clustering <- stats::hclust(dist_matrix)
  metabolites = order(metabolites)
  ordered_co_occurrence_matrix <- co_occurrence_matrix[metabolites, metabolites]
  diag(ordered_co_occurrence_matrix) = NA

  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(ordered_co_occurrence_matrix,
                                               name = "Co-occurrence",
                                               col = grDevices::colorRampPalette(c("white", "red"))(100),
                                               cluster_rows = FALSE,
                                               cluster_columns = FALSE,
                                               show_row_dend = FALSE,
                                               show_column_dend = FALSE,
                                               row_names_side = "left",  # Display row names on the left
                                               row_names_gp = grid::gpar(fontsize = 15),  # Increase row names font size
                                               column_names_gp = grid::gpar(fontsize = 15),  # Increase column names font size
                                               heatmap_legend_param = list(title = "Individuals flagged", legend_direction = "horizontal", title_position = "topcenter", title_gp = grid::gpar(fontsize = 15), labels_gp = grid::gpar(fontsize = 15)),  # Increase legend text size
                                               cell_fun = function(j, i, x, y, width, height, fill) {
                                                 if (grid::convertUnit(x, "native", valueOnly = TRUE)<=1-grid::convertUnit(y, "native", valueOnly = TRUE) | round(grid::convertUnit(x, "native", valueOnly = TRUE)+grid::convertUnit(y, "native", valueOnly = TRUE),2)==1) { # Only display the lower triangle
                                                   grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = fill, col = NA))
                                                   grid::grid.text(ifelse(is.na(ordered_co_occurrence_matrix[i, j])==FALSE,sprintf("%.0f", ordered_co_occurrence_matrix[i, j]),"NA"), x, y, gp = grid::gpar(fontsize = 12))
                                                 } else {
                                                   grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = "white", col = NA))
                                                 }
                                               }),
                       heatmap_legend_side = "top")

}
