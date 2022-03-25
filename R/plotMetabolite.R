`%>%` <- magrittr::`%>%`
#' @title
#' Initiate a data distribution visualization for the given metabolite.
#' @description
#' It shows a violin plot for different timepoints for the passed metabolites. Violin plot is a plot that shows a distribution of the
#' numeric list.
#'
#' @param data An object of class data.frame (or one that can be coerced to that class) containing data of all variables used in the plot
#' @param metabolite A character or a list of characters containing metabolites name for which the plot will be drawn. If a single metabolite is passed
#' then a plot for that metabolite is drawn. If a list of metabolite is passed then each metabolite will be shown in a different plot.
#' @param points As violin plot shows the distribution, observations along the end of tail are also shown by default. If you want to remove these points, then set points = FALSE.
#'
#' @return It returns plot for the metabolites passed.
#' @export
#'
#' @examples
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' plotMetabolite(data = metabol.data,metabolite = metabolites[1])
#' plotMetabolite(data = metabol.data,metabolite = metabolites[1:3],points = TRUE)
#'
plotMetabolite <- function(data, metabolite,points = FALSE){
  col_list = colnames(data)
  #Check if metabolite passed is present in data or not
  if(!is.character(metabolite) & any(sapply(metabolite, function(x){grepl(x,col_list)}))){
    stop("Metabolite passed is not present in the data")
  }

  if (length(metabolite)>1){
    for (item in metabolite) {
      plot_list = na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",item,".*",sep = ""))))
      new_data = reshape2::melt(data,measure.vars = plot_list,variable.name = "occurrence",value.name = "measurement")
      print(new_data %>% plotly::plot_ly(x=~occurrence,y=~measurement,split=~occurrence,type="violin",points = points) %>% plotly::layout(title = paste("Density of different timepoints for", item),
                                                                                                                    xaxis=list(title="Timepoints"),yaxis=list(title= paste(item, " Value"))))

    }
  }
  else{
    plot_list = na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolite,".*",sep = ""))))
    new_data = reshape2::melt(data,measure.vars = plot_list,variable.name = "occurrence",value.name = "measurement")
    new_data %>% plotly::plot_ly(x=~occurrence,y=~measurement,split=~occurrence,type="violin",points = points) %>% plotly::layout(title = paste("Density of different timepoints for", metabolite),
                                                                                                  xaxis=list(title="Timepoints"),yaxis=list(title= paste(metabolite, " Value")))
  }
}
