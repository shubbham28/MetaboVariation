#' @title
#' Extract the names of metabolites present in the data.
#' @description
#' Identify the unique metabolites present in the data. When metabolites are measured at different time points, the function returns a list of unique metabolites.
#'
#' @param list A list containing the names of columns that contain metabolite values across time points.
#' @param divider The character used to separate the metabolite name and time point name in the data. The default is "_".
#' @param start A binary value indicating if the metabolite name is followed by the time point or if the time point is followed by the metabolite.
#' For example, a column for metaboliteA at time point 1 named such as "MetaboliteA_1" should have start=TRUE, while a column name like "1_MetaboliteA" should have start=FALSE.
#'
#' @return A list containing the names of unique metabolites
#' @export
#'
#' @examples
#' # Load the simulated data.
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#'
#' # Extract the names of metabolites.
#' metabolites = get.metabolites(list = metabolite_list)
#' metabolites
get.metabolites <- function(list,divider = "_", start = TRUE){
  #Data check
  if(!(is.character(list) & is.character(divider))){
    stop("List and divider passed are not characters")
  }
  #If divider in the list or not
  if(!all(grepl(divider,list))){
    stop(paste("Following names does not contain the divider passed:",paste(list[!grepl(divider,list)],collapse = ", ")))
  }
  if (start == TRUE){
    sort(unique(sapply(list, function(x) {stringr::str_split(x,divider)[[1]][1]})))
  }
  else if(start == FALSE){
    sort(unique(sapply(list, function(x) {stringr::str_split(x,divider)[[1]][2]})))
  }
}
