#' @title
#' Extract metabolite names present in data
#' @description
#' get_metabolites() is a function that will extract the unique metabolites names from the columns present in data, when different timepoints of metabolites are listed as different columns in data.
#'
#' @param list A list containing names of columns that store metabolites values across timepoints.
#' @param divider A character value that separates the metabolite name and timepoint value. The default is "_"
#' @param start A binary value that tells the function if the metabolite name is followed by timepoint or the timepoint is followed by metabolite.
#' For eg: With divider "_", column names in format "MetaboliteA_1" should have start=TRUE, while column names in format "1_MetaboliteA" should have start=FALSE.
#'
#' @return A list of unique metabolite names
#' @export
#'
#' @examples
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' metabolites
get.metabolites <- function(list,divider = "_",start = TRUE){
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
