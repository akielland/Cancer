# importing and creating various dataframes and matrices


##############################################
## ALL DATA as they came original (I think) ##
##############################################

path <- "/Users/anders/Documents/MASTER/data/CORALLEEN_02.csv"
df01 <- read.table(path, header = TRUE, sep = ",", dec = ".")

###################################
## ALL DATA with added NODE DATA ##
###################################

# Adding node values to all data
path <- "/Users/anders/Documents/MASTER/data/model_predictors.tsv"
df_nodes <- read.delim(path)
colnames(df_nodes)[1] <- colnames(df01)[1]
df02 = full_join(df_nodes, df01, by = c(colnames(df01)[1]))


############################################
## Selecting data NODE DATA as predictors ##
############################################

df08 <- df07 |> filter(timepoint=="SUR") |>
  filter(TrialArmNeo=="Letro+Ribo") |>
  select(-1) |>
  select(-c(10:12))


#############
## 6 GENES ##
#############

df_6genes_with_output <- function(output){
  #Y_6genes <- output

  df_genes <- df01 |> select(c(UniqueID, CCND1, CCNE1, CDKN1A, ESR1, MYC, RB1, timepoint, TrialArmNeo)) |>
    filter(timepoint=="SCR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_output <- df01 |> select(c(UniqueID, output, timepoint, TrialArmNeo)) |>
    filter(timepoint=="SUR") |> 
    filter(TrialArmNeo=="Letro+Ribo")
  
  df_6genes <- full_join(df_genes, df_output, by = "UniqueID")
  #rename(df_6genes, c("Y") = c(output))
  colnames(df_6genes)[which(names(df_6genes) == output)] <- "Y"
  print(output)
  df_6genes <- na.omit(df_6genes)
  return(df_6genes)
}
  
Proliferation_6genes <- as.data.frame(df_6genes_with_output("ProliferationScore"))




####################################################################
## Prevoius code:
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
# Y from tabel with the node values
Y <- select(df08, "ProliferationScore")

#########################
##   ##
#########################


