# MetaboVariation

Before we start the MetaboVariation walkthrough, the user should have working R software enviroment installed on their machine. The MetaboVaraiton package has following dependencies and will also be installed along with the package if not already installed on the machine.

* brms
*  circlize
*  ComplexHeatmap
*  doParallel
*  dplyr
*  foreach
*  future
*  grid
*  magrittr
*  parallel
*  plotly
*  reshape2
*  rstan
*  scales
*  stringr
*  readxl
*  tidyr
*  ggplot2

To install this package, start by installing the devtools package. The best way to do this is from CRAN, by typing:

**install.packages("devtools")**

After successful installation, use the following code to install the package on your system:

**library(devtools)**

**install_github("shubbham28/MetaboVariation")**
