# importing and creating various dataframes and matrices


# df01: 6 genes; timepoint SCR and SUR
df01 <- read.delim2('/Users/anders/Documents/MASTER/Cancer/R_codeP01/data/processed/SCR_SUR_6genes_noNAN(ANO).txt')
df01 <- select(df01, 1:9)

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
