


dat <- list(ROR_prolif_771genes=ROR_prolif_771genes, char_list=char_list, synergistic=synergistic)
save(dat, file = "dat.RData")


a=names(result$beta_main)
lapply(char_list, function(xx) sum(a %in% as.vector(xx)))


sum(a %in% as.vector(char_list$prolif_))

as.vector(char_list$prolif_)

a=names(a)
a
# [1] "ALDOA"    "AURKB"    "AXIN2"    "BAIAP2L1" "BBOX1"    "BMP7"     "BMPR1B"  
# [8] "CACNA1H"  "CD44"     "COL9A3"   "DDX39A"   "DKK1"     "EIF3B"    "EP300"   
# [15] "EYA1"     "FGF12"    "FOXC1"    "FZD10"    "HBB"      "HDAC5"    "LEF1"    
# [22] "LIF"      "MLH1"     "MMP3"     "PLAT"     "PMS2"     "PRLR"     "PROM1"   
# [29] "RBL2"     "RBX1"     "RELN"     "RORB"     "SHMT2"    "SP1"      "TLE3"    
# [36] "TUBA4A"   "UBE2T"   

sum(a %in% as.vector(char_list$prolif_))
[1] 5

sum(a %in% as.vector(char_list$ER_sing_))
[1] 4





