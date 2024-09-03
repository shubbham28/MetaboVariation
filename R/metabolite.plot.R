`%>%` <- magrittr::`%>%`
#' @title
#' Plot the distribution of metabolite levels.
#' @description
#' A function that produces violin plots of the distributions of metabolites for all individuals across all time points.
#'
#' @param data A data frame containing data of all variables to be used in the BGLM. Refer to \code{\link{metabol.data}} for the required structure and content of the data frame.
#' @param metabolite A string or a list of strings containing the names of the metabolites to be plotted.
#' @param title The text for the title of the plot.
#'
#' @export
#'
#' @examples
#' # Load the simulated data and extract the metabolites names.
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#'
#' # Plots the distribution of a single metabolite or multiple metabolites.
#' metabolite.plot(data = metabol.data,metabolite = metabolites[1])
#' metabolite.plot(data = metabol.data,metabolite = metabolites[1:3])
#'
metabolite.plot <- function(data, metabolite,title=''){
  col_list = colnames(data)
  #Check if metabolite passed is present in data or not
  if(!is.character(metabolite) & any(sapply(metabolite, function(x){grepl(x,col_list)}))){
    stop("Metabolite passed is not present in the data")
  }

  if (length(metabolite)>1){
    for (item in metabolite) {
      plot_list = stats::na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",item,".*",sep = ""))))
      new_data = reshape2::melt(data,measure.vars = plot_list,variable.name = "occurrence",value.name = "measurement")
      levels(new_data$occurrence) = paste0("Timepoint ",1:length(levels(new_data$occurrence)))
      print(new_data %>% plotly::plot_ly(x=~occurrence,y=~measurement,split=~occurrence,type="violin",points = FALSE) %>% plotly::layout(title = title,
                                                                                                                                         xaxis=list(title="Timepoints"),yaxis=list(title= paste("Value of ",item))))

    }
  }
  else{
    plot_list = stats::na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolite,".*",sep = ""))))
    new_data = reshape2::melt(data,measure.vars = plot_list,variable.name = "occurrence",value.name = "measurement")
    levels(new_data$occurrence) = paste0("Timepoint ",1:length(levels(new_data$occurrence)))
    new_data %>% plotly::plot_ly(x=~occurrence,y=~measurement,split=~occurrence,type="violin",points = FALSE) %>% plotly::layout(title = title,
                                                                                                                                 xaxis=list(title="Timepoints"),yaxis=list(title= paste("Value of ",metabolite)))
  }
}
