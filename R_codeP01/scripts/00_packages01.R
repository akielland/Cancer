library(MASS)

library(tidyverse)
library(readxl)
library(tidymodels)
library(dplyr)
library(ggplot2)

library(reshape2)

library(gridExtra)

library(glmnet)
library(boot)
library(caret)
library(mboost)

library(sparsgel)



install.packages("grpreg")
install.packages("sparsegl")

install.packages("gridExtra")


# Install and load the sparsgel package
if (!requireNamespace("sparsgel", quietly = TRUE)) {
  install.packages("sparsgel")
}



install.packages("caret", dependencies = c("Depends", "Suggests"))
# https://topepo.github.io/caret/


install.packages("ModelMetrics")

# install.packages("tidyverse")

# On the pop-up, "Do you want to install from source the packages which need compilation?", 
# [No] will ensure that all packages get updated, but not necessarily to their latest versions. 
# [Yes] should update everything to its latest version, but only if you installed the latest Rtools. 
# If Rtools is not up to date, then among packages that have not yet been compiled to binaries, [Yes] 
# will successfully update some (or none) and will fail on some (or none). 
# If any fail, then update.packages() can be run again, selecting [No] to get the latest versions 
# available without updating Rtools.

class(data_in$Y)
data1[, "Y"]
