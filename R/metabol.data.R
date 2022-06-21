#'  Metabolite data of individuals
#'
#' A simulated dataset containing metabolite levels of 3 metabolites for 164 individuals measured at four different time points. The covariates sex, age and BMI are also available for these individuals.
#'
#' @format A data frame with 164 rows and 16 columns. Missing values for the individual for that time point is showed as NA.
#' \describe{
#'   \item{Age}{The age of the individuals.}
#'   \item{BMI}{The BMI of the individuals.}
#'   \item{SexM.1F.2}{The gender of the individuals. It is a factor where 1 denotes "Male" while 2 denotes "Female".}
#'   \item{metabolA_1}{The value of metabolite "A" for the first timepoint.}
#'   \item{metabolA_2}{The value of metabolite "A" for the second timepoint.}
#'   \item{metabolA_3}{The value of metabolite "A" for the third timepoint.}
#'   \item{metabolA_4}{The value of metabolite "A" for the fourth timepoint.}
#'   \item{metabolB_1}{The value of metabolite "B" for the first timepoint.}
#'   \item{metabolB_2}{The value of metabolite "B" for the second timepoint.}
#'   \item{metabolB_3}{The value of metabolite "B" for the third timepoint.}
#'   \item{metabolB_4}{The value of metabolite "B" for the fourth timepoint.}
#'   \item{metabolC_1}{The value of metabolite "C" for the first timepoint.}
#'   \item{metabolC_2}{The value of metabolite "C" for the second timepoint.}
#'   \item{metabolC_3}{The value of metabolite "C" for the third timepoint.}
#'   \item{metabolC_4}{The value of metabolite "C" for the fourth timepoint.}
#'   \item{Individual_id}{The id of the individuals.}
#' }
#' @usage data(metabol.data)
"metabol.data"
