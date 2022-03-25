#'  Metabolites data of individuals
#'
#' A random dataset that contains metabolites levels of 3 metabolites for different individuals taken on different timepoints.
#'
#' @format A data frame with 164 rows and 16 columns. NA in the data means that the metabolite value is missing for the individual for that timepoint.
#' \describe{
#'   \item{Age}{This contains the age of the individual}
#'   \item{BMI}{This contains the BMI value of the individual}
#'   \item{SexM.1F.2}{Gender of the individual is recorded here. It's a factor column where 1 denotes "Male" while 2 denotes "Female"}
#'   \item{metabolA_1}{Contains the value of "A" metabolite for the first timepoint.}
#'   \item{metabolA_2}{Contains the value of "A" metabolite for the second timepoint.}
#'   \item{metabolA_3}{Contains the value of "A" metabolite for the third timepoint.}
#'   \item{metabolA_4}{Contains the value of "A" metabolite for the fourth timepoint.}
#'   \item{metabolB_1}{Contains the value of "B" metabolite for the first timepoint.}
#'   \item{metabolB_2}{Contains the value of "B" metabolite for the second timepoint.}
#'   \item{metabolB_3}{Contains the value of "B" metabolite for the third timepoint.}
#'   \item{metabolB_4}{Contains the value of "B" metabolite for the fourth timepoint.}
#'   \item{metabolC_1}{Contains the value of "C" metabolite for the first timepoint.}
#'   \item{metabolC_2}{Contains the value of "C" metabolite for the second timepoint.}
#'   \item{metabolC_3}{Contains the value of "C" metabolite for the third timepoint.}
#'   \item{metabolC_4}{Contains the value of "C" metabolite for the fourth timepoint.}
#'   \item{Individual_id}{Contains the individual id for the individual that provided the measurements for various metabolites.}
#' }
#' @usage data(metabol.data)
"metabol.data"
