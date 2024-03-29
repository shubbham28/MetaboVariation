`%>%` <- magrittr::`%>%`
#' @title
#' Plots the distribution of metabolite levels.
#' @description
#' A function that plots (as a violin plot) the distributions of metabolites for all individuals across all time points.
#'
#' @param data A data frame containing data on all individuals
#' @param metabolite A string or a list of strings containing the names of the metabolites to be plotted.
#' @param main The text for the main title.
#'
#' @return A violin plot for the specified metabolites.
#' @export
#'
#' @examples
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' metabolite.plot(data = metabol.data,metabolite = metabolites[1])
#' metabolite.plot(data = metabol.data,metabolite = metabolites[1:3])
#'
metabolite.plot <- function(data, metabolite,main=''){
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
      print(new_data %>% plotly::plot_ly(x=~occurrence,y=~measurement,split=~occurrence,type="violin",points = FALSE) %>% plotly::layout(title = main,
                                                                                                                                         xaxis=list(title="Timepoints"),yaxis=list(title= paste("Value of ",item))))

    }
  }
  else{
    plot_list = stats::na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolite,".*",sep = ""))))
    new_data = reshape2::melt(data,measure.vars = plot_list,variable.name = "occurrence",value.name = "measurement")
    levels(new_data$occurrence) = paste0("Timepoint ",1:length(levels(new_data$occurrence)))
    new_data %>% plotly::plot_ly(x=~occurrence,y=~measurement,split=~occurrence,type="violin",points = FALSE) %>% plotly::layout(title = main,
                                                                                                                                 xaxis=list(title="Timepoints"),yaxis=list(title= paste("Value of ",metabolite)))
  }
}
