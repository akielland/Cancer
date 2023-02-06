###############################################################
## Ensemble: create data structures from the signature genes ##
###############################################################

# old structures of base learners
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm05 <- as.formula(paste("Y", paste(genes, collapse="+"), sep=" ~ "))

## formulas from signatuer genes

signature <- function(words){
  split_words <- strsplit(words, "\t")    # if also space between words use;  "\t| ". 
  # The do.call function is used to combine the split words into a single character matrix
  # Convert each split word into a character vector and combine into a single character matrix
  character_matrix <- do.call(rbind, split_words)
  # create the formula
  fm_signature <- as.formula(paste("Y", paste(character_matrix, collapse="+"), sep=" ~ ")) 
  return(fm_signature)
}

## Immune Infiltration
words <- c("CD8A	HLA.DRB1	PSMB10	CCL5	APOE	CD27	CD274	CD276	CHIT1	CLEC5A	CMKLR1	COLEC12	CXCL5	CXCL9	CXCR6	CYBB	GZMA	GZMB	GZMH	GZMM	HLA.E	IDO1	LAG3	MARCO	MSR1	NKG7	PDCD1	PDCD1LG2	SIGIRR	SPN	STAT1	TIGIT	ZNF205")
fm_immune_inf <- signature(words)
fm_immune_inf

## Proliferation
words <- c("ITGB3	MMP9	CCNE1	EDN1	HIF1A	KDR	TYMP	ANLN	APOD	BMP7	CD24	CDKN2A	CXCL12	CXCL13	FGFR3	IBSP	MAPK1	MAPK3	THBS4	ACVR1C	ACVRL1	AGT	AREG	ASPN	ATM	AURKA	AURKB	BAD	BMP4	BMP5	BMPR1A	BMPR1B	BMPR2	CBLC	CCNA1	CCNA2	CCNB1	CCND1	CCND2	CCNE2	CCR2	CDC14A	CDC14B	CDC20	CDC25A	CDC25B	CDC6	CDC7	CDK4	CDK6	CDKN1A	CDKN1B	CDKN1C	CDKN2B	CDKN2C	CDKN2D	CENPF	CEP55	CETN2	CHAD	CHEK2	CREBBP	E2F1	E2F5	EGF	EGFR	EGLN2	ENO1	EP300	ERBB2	EREG	EXO1	F3	FGF2	FGF9	FGFR2	FGFR4	FHL1	GADD45A	GADD45B	GADD45G	GATA4	GDF15	GDF5	GREM1	GSK3B	HDAC1	HDAC2	HGF	ID1	IGF1	IL1B	IL6	INHBA	JAG1	KIF2C	LIFR	MCM2	MDM2	MELK	MET	MKI67	MTOR	MYC	NDP	NPR1	OGN	PCNA	PDGFRB	PIK3CA	PIK3R2	PKMYT1	PPP2R1A	PRKCB	PRKDC	PTEN	PTGS2	RB1	RBX1	RPS6KB2	RRM2	SFN	SFRP2	SKP1	SKP2	SMAD1	SMAD3	SMAD4	SMAD5	SMC1B	SMURF2	SOX17	TFDP1	TGFB1	TGFB2	TGFB3	TGFBR2	TNN	TP53	UBE2C	VEGFA	VEGFD	WEE1	ZFYVE9")
fm_prolif <- signature(words)
fm_prolif

## Angiogenesis	
words <- c("ANGPT1	CD44	CDH5	ITGB3	MMP9	TEK	VCAN	BCL6B	CCNE1	CLEC14A	CXorf36	EDN1	EMCN	EPAS1	FAM124B	FGF18	FSTL3	HIF1A	IKZF3	ITGAV	
           KDR	MMRN2	MYCT1	PALMD	PDGFB	RNASE2	ROBO4	SERPINB5	SERPINH1	SHE	STC1	TIE1	TNFAIP6	TYMP")
fm_angiogenesis <- signature(words)
fm_angiogenesis

## Antigen Presentation	
words <- c("THBS1	CD1E	CD8A	HLA-A	HLA-B	HLA-C	HLA-DMA	HLA-DMB	HLA-DOB	HLA-DPA1	HLA-DPB1	HLA-DQA1	HLA-DQB1	HLA-DRA	HLA-DRB1	PSMB10	PSMB7	PSMB9	TAP1	TAP2	TAPBP")
fm_antigen_present <- signature(words)
fm_antigen_present

## ER Signaling	
words <- c("ADCY9	ADD1	ANXA9	BORCS7	CDCA8	DDX39A	DNAJC12	EIF3B	ELOVL2	ESR1	HEMK1	IFT140	ITPR1	MAPT	NAT1	PFDN2	PGR	PTGER3	SCUBE2	SERBP1	SHMT2	SYTL4	TBC1D9	TCEAL1	TFF1	TLE3	WDR77")
fm_ER_signaling <- signature(words)
fm_ER_signaling





