#' @title
#' Extract names of metabolites present in data
#' @description
#' Identify the unique metabolites present in the data. When metabolites are present at different time points, the function returns a list of unique metabolites.
#'
#' @param list A list containing the name of columns that contain metabolite values across time points.
#' @param divider A character separating the metabolite name and time point name. The default is "_"
#' @param start A binary value indicating if the metabolite name is followed by the time point or if the time point is followed by metabolite.
#' For example, column name such as "MetaboliteA_1" should have start=TRUE, while if column name like "1_MetaboliteA" should have start=FALSE.
#'
#' @return A list containing the names of unique metabolites
#' @export
#'
#' @examples
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = getmetabolites(list = metabolite_list)
#' metabolites
getmetabolites <- function(list,divider = "_",start = TRUE){
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
