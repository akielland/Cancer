##############################################################
## Importing and creating various dataframes and matrices 
##
## Here: focusing on testing
# - ROC
# - 2.trail


library(tidyverse)
library(readxl)


########################################################
## ALL DATA from the second trail (a test set for us) ##
########################################################

path <- "~/Documents/Master/data/COR_and_validation_data_with_model_prediction/cdk_cohort_with_model_prediction.tsv"
df04 <- read.delim(path)


# extracting genes
ext.X <- df04 |> 
  select(42:812)

dim(ext.X)
ext.X <- as.matrix(ext.X)

# extracting progress free survival (PFS)
ext.PFS_months <- df04 |> select(PFS_months)
ext.PFS_status <- df04 |> select(PFS_status)




########################################
## 771 GENES + ROR categorical output ##
########################################

# select all gene names from original data table structure
genes <- colnames(df01 |> select(41:811))
length(genes)


df_genes_ROR_classification <- function(df, predictors, output){
  features <- c("UniqueID", "timepoint", "TrialArmNeo", predictors)
  
  df_genes <- df |> select(all_of(features)) |>
    filter(timepoint=="SCR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_SUR <- df |> select(c(UniqueID, timepoint, TrialArmNeo, output, ROR_P_GroupSubtypeProliferation)) |>
    filter(timepoint=="SUR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_out <- full_join(df_genes, df_SUR, by = "UniqueID") |> 
    select(-UniqueID)
  
  colnames(df_out)[which(names(df_out) == output)] <- "Y"
  print(output)
  df_out <- na.omit(df_out) # remove rows with missing values
  return(df_out)
}

df.ROR <- df_genes_ROR_classification(df01, genes, "ROR_P_Subtype_Proliferation") 

# decode ROR "high" => 1 and "med"/"low" to 0 
df.ROR <- df.ROR |>  
  mutate(ROR_P_GroupSubtypeProliferation = ifelse(ROR_P_GroupSubtypeProliferation == "high", 1, 0))
# some testing of for the code above
test <- df.ROR |>  
  mutate(ROR_P_Group = ifelse(ROR_P_GroupSubtypeProliferation %in% c("high"), 1, 0))
test <- df.ROR |>  
  mutate(ROR_P_Group = ifelse(ROR_P_GroupSubtypeProliferation == "high", 1, 0))
test |> select(ROR_P_Group, ROR_P_GroupSubtypeProliferation)



# substitute in the data frame (make new dataframe with scaled values of pred)
df.all_genes.pred_mech.scaled <- df.all_genes.pred_mech |> 
  mutate(model_prediction = mech_pred_scaled)

# check if data looks the same
hist(df.all_genes.pred_mech.scaled$model_prediction)


# Matrices for lasso
X <- as.matrix(df.all_genes.pred_mech.scaled |> 
                 select(genes))
pred_mech <- as.matrix(select(df.all_genes.pred_mech.scaled, model_prediction))
Y <- as.matrix(select(df.all_genes.pred_mech.scaled, Y))












