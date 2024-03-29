% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.metabolites.R
\name{get.metabolites}
\alias{get.metabolites}
\title{Extract names of metabolites present in data}
\usage{
get.metabolites(list, divider = "_", start = TRUE)
}
\arguments{
\item{list}{A list containing the name of columns that contain metabolite values across time points.}

\item{divider}{The character used to separate the metabolite name and time point name in the data. The default is "_"}

\item{start}{A binary value indicating if the metabolite name is followed by the time point or if the time point is followed by the metabolite.
For example, a column name such as "MetaboliteA_1" should have start=TRUE, while a column name like "1_MetaboliteA" should have start=FALSE.}
}
\value{
A list containing the names of unique metabolites
}
\description{
Identify the unique metabolites present in the data. When metabolites are present at different time points, the function returns a list of unique metabolites.
}
\examples{
data(metabol.data)
metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
metabolites = get.metabolites(list = metabolite_list)
metabolites
}
