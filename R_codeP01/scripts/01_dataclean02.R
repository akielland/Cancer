# Data for machine predictions from nodes

# importing and creating various dataframes and matrices

#############
## NODES ##
#############

path_nodes <- "/Users/anders/Documents/MASTER/data/For_NODES.tsv"

df05 <- read.delim(path_nodes)
df05 <- select(df05, -1)
df05 <- select(df05, -10)

# df02: df01 but done wider by removimin the timpoint colum and make uniqe names for features at each timepiont
# ex scr_CCNED1 and sur_CCNED1
df01_scr <- df01 |> filter(timepoint=="SCR")
colnames(df01_scr) <- paste("scr", colnames(df01_scr), sep = "_")
df01_sur <- df01 |> filter(timepoint=="SUR")
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



#########################
##   ##
#########################


