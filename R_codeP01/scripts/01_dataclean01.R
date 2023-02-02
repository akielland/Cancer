############################################################
## importing and creating various dataframes and matrices ##
############################################################

# REMEMBER: Standardizing features can be relevant
# it is extremely important that all PREDICTOR variables are on the same scale. 
# These regularization methods are sensitive to the size of coefficients by design, 
# so coefficients on different scales are enormously problematic. 
# Based on how these data were simulated, we know that our predictors are already on the same scale, 
# but as good practice, we will center and scale them anyways.
# Use scale()
pred <- scale(pred)
X_Z <- scale(X)

library(tidyverse)
library(readxl)

##############################################
## ALL DATA as they came original (I think) ##
##############################################
path <- "/Users/anders/Documents/MASTER/data/CORALLEEN_02.csv"
df01 <- read.table(path, header = TRUE, sep = ",", dec = ".")

##########################################################
## ALL DATA with added NODE DATA from MECHANISTIC MODEL ##
##########################################################
# Adding node values to all data
path <- "/Users/anders/Documents/MASTER/data/model_predictors.tsv"
df_nodes <- read.delim(path)
colnames(df_nodes)[1] <- colnames(df01)[1]
df02 = full_join(df_nodes, df01, by = c(colnames(df01)[1]))


#############################################################
## ALL data with added prediction of the Mechanistic model ##
#############################################################
path <- "/Users/anders/Documents/Master/data/COR_and_validation_data_with_model_prediction/CORALLEEN_data_with_model_prediction.xlsx"
# df_03 <- read_excel(path, range = cell_cols("A:V"))
df03 <- read_excel(path) |> 
  rename(TrialArmNeo = "Trial Arm Neo")

############################################################
## ALL data with added residuals of the Mechanistic model ##
############################################################

# prediction of proliferation.score
# residuals = ??
path <- "/Users/anders/Documents/Master/data/COR_and_validation_data_with_model_prediction/CORALLEEN_data_with_model_prediction.xlsx"
# df_03r <- read_excel(path, range = cell_cols("A:V")) |> mutate(residuals = Proliferation.Score - model_prediction)
df03r <- read_excel(path) |>
  mutate(residuals = Proliferation.Score - model_prediction) |> 
  rename(TrialArmNeo = "Trial Arm Neo")

#############
## 6 GENES ##
#############

df_6genes_with_output <- function(output){
  print(output)
  df_genes <- df01 |> select(UniqueID, CCND1, CCNE1, CDKN1A, ESR1, MYC, RB1, timepoint, TrialArmNeo) |>
    filter(timepoint=="SCR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_output <- df01 |> select(c(UniqueID, output, timepoint, TrialArmNeo)) |>
    filter(timepoint=="SUR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_6genes <- full_join(df_output, df_genes, by = "UniqueID")
  
  #rename(df_6genes, c("Y") = c(output))
  colnames(df_6genes)[which(names(df_6genes) == output)] <- "Y"
  
  df_6genes <- df_6genes |> select(-c(UniqueID, timepoint.x, TrialArmNeo.x, timepoint.y, TrialArmNeo.y))
  df_6genes <- na.omit(df_6genes) # remove rows with NA
  return(df_6genes)
}

prolif_6genes <- df_6genes_with_output("ProliferationScore")
RORprolif_6genes <- df_6genes_with_output("ROR_P_Subtype_Proliferation")


# Matrices for lasso and XGboost (not needed anymore/done in functions)
X <- as.matrix(Proliferation_6genes |> select(CCND1, CCNE1, CDKN1A, ESR1, MYC, RB1))
Y <- as.matrix(select(Proliferation_6genes, Y))

###############
## 771 GENES ##
###############
# This code snip is more generic and can be used instead off others i think

# select ALL gene names from original data table structure
all_genes <- colnames(df01 |> select(41:811))

# just for testing with fewer features
test_genes <- colnames(df01 |> select(41:46))

length(genes)

df_genes_with_output <- function(df, predictors, output){
  print(output)
  print(length(genes))
  features <- c("UniqueID", "timepoint", "TrialArmNeo", predictors)
  
  df_genes <- df |> select(all_of(features)) |>
    filter(timepoint=="SCR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_output <- df |> select(c(UniqueID, output, timepoint, TrialArmNeo)) |>
    filter(timepoint=="SUR") |> 
    filter(TrialArmNeo=="Letro+Ribo")

  df_out <- full_join(df_output, df_genes, by = "UniqueID")
  
  colnames(df_out)[which(names(df_out) == output)] <- "Y"
  
  df_out <- df_out |> select(-c(UniqueID, timepoint.x, TrialArmNeo.x, timepoint.y, TrialArmNeo.y)) # remove unnecessary columns
  df_out <- na.omit(df_out) # remove rows with missing values
  return(df_out)
}

prolif_771genes <- df_genes_with_output(df01, all_genes, "ProliferationScore")
RORprolif_771genes <- df_genes_with_output(df01, all_genes, "ROR_P_Subtype_Proliferation")

head(prolif_771genes[,1:5])
lastcol <- ncol(prolif_771genes)
head(prolif_771genes[, (lastcol-5): lastcol])

# dataframe for caret
dfA03 <- Proliferation_ALLgenes |> 
  select(- c(timepoint.x, TrialArmNeo.x, timepoint.y, TrialArmNeo.y))

# Matrices for glmnet (not in use anymore)
X <- as.matrix(Proliferation_ALLgenes |> select(all_of(genes)))
Y <- as.matrix(select(Proliferation_ALLgenes, Y))


#################
## NODES DATA  ##
#################
# Dataframe with values of the nodes by the mechanistic model

df_Nodes_with_output <- function(output){
  df_features <- df02 |>
    select(UniqueID, cyclinD1, cyclinD1Palbo, p21, cyclinD1p21, cMyc, cyclinEp21, Rb1, ppRb1, timepoint, TrialArmNeo) |> 
    filter(timepoint=="SCR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_output <- df02 |> select(UniqueID, output, timepoint, TrialArmNeo) |>
    filter(timepoint=="SUR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_Nodes_Y <- full_join(df_output, df_features, by = "UniqueID")
  
  colnames(df_Nodes_Y)[which(names(df_Nodes_Y) == output)] <- "Y"
  print(output)
  
  df_Nodes_Y <- df_Nodes_Y |> select(-c(UniqueID, timepoint.x, TrialArmNeo.x, timepoint.y, TrialArmNeo.y)) # remove unnecessary columns
  df_Nodes_Y <- na.omit(df_Nodes_Y) # remove rows containing NA
  return(df_Nodes_Y)
}

prolif_nodes <- df_Nodes_with_output("ProliferationScore")
RORprolif_nodes <- df_Nodes_with_output("ROR_P_Subtype_Proliferation")

# Matrices for lasso (not used for newer code)
X <- as.matrix(Nodes_Proliferation |> 
                 select(cyclinD1, cyclinD1Palbo, p21, cyclinD1p21, cMyc, cyclinEp21, Rb1, ppRb1))
Y <- as.matrix(select(Nodes_Proliferation, Y))

#################################
## 771 GENES + MECH prediction ##
#################################
# select all gene names from original data table structure
genes <- colnames(df03 |> select(41:811))
length(genes)

df_genes_with_mech.pred_and_output <- function(df, predictors, output){
  features <- c("UniqueID", "timepoint", "TrialArmNeo", predictors)
  
  df_genes <- df |> select(all_of(features)) |>
    filter(timepoint=="SCR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_SUR <- df |> select(c(UniqueID, timepoint, TrialArmNeo, output, model_prediction)) |>
    filter(timepoint=="SUR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_out <- full_join(df_SUR, df_genes, by = "UniqueID") |> 
    select(-UniqueID)
  
  df_out <- df_out |> 
    select(- c(timepoint.x, TrialArmNeo.x, timepoint.y, TrialArmNeo.y))  # remove unnecessary columns
  
  colnames(df_out)[which(names(df_out) == output)] <- "Y"
  print(output)
  df_out <- na.omit(df_out) # remove rows with missing values
  return(df_out)
}

df_771genes_mech_pred_prolif <- df_genes_with_mech.pred_and_output(df03, genes, "Proliferation.Score") 

## scale mech model prediction to scale of proliferation score by a linear model
# fit model
fit_ <- lm(df_771genes_mech_pred_prolif$Y ~ df_771genes_mech_pred_prolif$model_prediction)
summary(fit_)
plot(df_771genes_mech_pred_prolif$model_prediction, df_771genes_mech_pred_prolif$Y)
abline(fit_ , col = "blue")

# scale value of mechanistic model by the linear model
mech_pred_scaled <- predict(fit_, newdata =  as.data.frame(df_771genes_mech_pred_prolif$model_prediction))
cor(df_771genes_mech_pred_prolif$model_prediction, mech_pred_scaled)
hist(df_771genes_mech_pred_prolif$model_prediction)
cor(df_771genes_mech_pred_prolif$Y, df_771genes_mech_pred_prolif$model_prediction, method="pearson")
cor(df_771genes_mech_pred_prolif$Y, df_771genes_mech_pred_prolif$model_prediction, method="spearman")
plot(mech_pred_scaled, df_771genes_mech_pred_prolif$Y)

# substitute in the data frame (make new data frame with scaled values of prediction)
df_771genes_mech_pred.scaled_prolif <- df_771genes_mech_pred_prolif |> 
  mutate(model_prediction = mech_pred_scaled)

# check if data looks the same
hist(df_771genes_mech_pred.scaled_prolif$model_prediction)


# Matrices for lasso (not necessary for newer code)
X <- as.matrix(df_771genes_mech_pred.scaled_prolif |> select(genes))
pred_mech <- as.matrix(select(df_771genes_mech_pred.scaled_prolif, model_prediction))
Y <- as.matrix(select(df_771genes_mech_pred.scaled_prolif, Y))


#####################
## OLD code:
#####################

# df01: 6 genes; timepoint SCR and SUR
path_6genes <- "/Users/anders/Documents/MASTER/data/SCR_SUR_6genes_noNAN(ANO).txt"

df01 <- read.delim2(path_6genes)
df01 <- select(df01, 1:9)

# df02: df01 but done wider by removimin the timpoint colum and make uniqe names for features at each timepiont
# ex scr_CCNED1 and sur_CCNED1
df01_scr <- df01 |> 
  filter(timepoint=="SCR")
colnames(df01_scr) <- paste("scr", colnames(df01_scr), sep = "_")
df01_sur <- df01 |> 
  filter(timepoint=="SUR")
colnames(df01_sur) <- paste("sur", colnames(df01_sur), sep = "_")

df02 <- cbind(df01_scr, df01_sur) |> 
  select(-ends_with("timepoint")) |> 
  select(-"sur_TrialArmNeo") |> 
  rename("TrialArmNeo" = "scr_TrialArmNeo")

# divided df02 into 2 new df with each trial arm
df02_LR <- filter(df02, TrialArmNeo=="Letro+Ribo")
df02_AC <- filter(df02, TrialArmNeo=="AC+pacli")


###############
## ALL GENES ##
###############

path_ALL <- "/Users/anders/Documents/MASTER/data/CORALLEEN_minus_3.txt"

df03 <- read.delim2(path_ALL, check.names = F)

df04 <- df03 |> 
  select(-c(6:13)) |> 
  select(-c(7:15)) |> 
  select(-c(8:16)) |> 
  filter(TrialArmNeo=="Letro+Ribo")
  
df04_scr_genes <- df04 |>   
  filter(timepoint=="SCR") |> 
  select(-c(1:7))

df04_sur_prolif <- df04 |> 
  filter(timepoint=="SUR") |> 
  select("Proliferation.Score", "Ki67")

df04 <- cbind(df04_sur_prolif, df04_scr_genes)

#########################
## X and Y's for lasso ##
#########################

# X and Y from the 6 gene data
X <- dplyr::select(df02_LR,  scr_CCND1, scr_CCNE1, scr_CDKN1A, scr_ESR1, scr_MYC, scr_RB1)
Y <- select(df02_LR, sur_ProliferationScore)

# X from the full data set
X <- select(df04, -"Proliferation.Score") |> 
  select(-"Ki67")

# Y's from the full data set
Y <- select(df04, "Proliferation.Score")
Y <- select(df04, "Ki67")
Y[is.na(Y)] = 0  # change NA too 0

# X from the node values
X <- dplyr::select(df08, cyclinD1, cyclinD1Palbo, p21, cyclinD1p21, cMyc, cyclinEp21, Rb1, ppRb1, proliferation)
X <- dplyr::select(df08, cyclinD1, cyclinD1Palbo, p21, cyclinD1p21, cMyc, cyclinEp21, Rb1, ppRb1)
# Y from table with the node values
Y <- select(df08, "ProliferationScore")



#########################
##   ##
#########################


