#'  Metabolite data of individuals
#'
#' A simulated dataset containing metabolite levels of 3 metabolites for 150 individuals measured at four different time points. The covariates sex, age and BMI are also available for these individuals.
#'
#' @format A data frame with 150 rows and 16 columns. Missing values for the individual for that time point is showed as NA.
#' \describe{
#'   \item{Age}{The age of the individuals.}
#'   \item{BMI}{The BMI of the individuals.}
#'   \item{SexM.1F.2}{The gender of the individuals. It is a factor where 1 denotes "Male" while 2 denotes "Female".}
#'   \item{metaboliteA_1}{The value of metabolite "A" for the first timepoint.}
#'   \item{metaboliteA_2}{The value of metabolite "A" for the second timepoint.}
#'   \item{metaboliteA_3}{The value of metabolite "A" for the third timepoint.}
#'   \item{metaboliteA_4}{The value of metabolite "A" for the fourth timepoint.}
#'   \item{metaboliteB_1}{The value of metabolite "B" for the first timepoint.}
#'   \item{metaboliteB_2}{The value of metabolite "B" for the second timepoint.}
#'   \item{metaboliteB_3}{The value of metabolite "B" for the third timepoint.}
#'   \item{metaboliteB_4}{The value of metabolite "B" for the fourth timepoint.}
#'   \item{metaboliteC_1}{The value of metabolite "C" for the first timepoint.}
#'   \item{metaboliteC_2}{The value of metabolite "C" for the second timepoint.}
#'   \item{metaboliteC_3}{The value of metabolite "C" for the third timepoint.}
#'   \item{metaboliteC_4}{The value of metabolite "C" for the fourth timepoint.}
#'   \item{Individual_id}{The id of the individuals.}
#' }
#' @usage data(metabol.data)
"metabol.data"
