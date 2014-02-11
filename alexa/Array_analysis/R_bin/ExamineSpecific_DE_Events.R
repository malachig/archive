#Written by Malachi Griffith
#Load neccessary libraries
library(multtest)
library(nortest)

#1.) Import DATA
#Load exon expression values for each probe
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)



exon_exp_data = read.table(file="Standard_probe_data.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(4))
length(exon_exp_data[,1])

x = which((exon_exp_data[,"Probe_Type"] == "Exon") | (exon_exp_data[,"Probe_Type"] == "Exon-Exon" & exon_exp_data[,"Exons_Skipped"] == 0))
canonical_exp_data = exon_exp_data[x,]
length(canonical_exp_data[,1])


######################################################################################
#ALEXA Gene: 5015 - AKAP18
gene_id = 5015
de_isoform_probesets = c(2388201,664004,2388202,664006,2388203)
length(de_isoform_probesets)
de_isoform_data = canonical_exp_data[which(canonical_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])



######################################################################################
#ALEXA Gene: 21000 - C12orf63
gene_id = 21000
de_isoform_probesets = c(1819129,1819127,2508111,2508112,1819085,2508110,2508109,2508103,1819124,2508106,1819109,1819120,2508104,2508105,1819102,2508107,1819115,1819094,2508108)
length(de_isoform_probesets)
de_isoform_data = canonical_exp_data[which(canonical_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])

######################################################################################
#ALEXA Gene: 3489 - C5
gene_id = 3489
de_isoform_probesets = c(2370135,2370134,486942)
length(de_isoform_probesets)
de_isoform_data = canonical_exp_data[which(canonical_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 2569 - CDC25B
gene_id = 2569
de_isoform_probesets = c(2360015)
length(de_isoform_probesets)
de_isoform_data = canonical_exp_data[which(canonical_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 5830 - COL21A1
gene_id = 5830
de_isoform_probesets = c(2095113,2095112,2095110)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 1645 - DIS3
gene_id = 1645
de_isoform_probesets = c(282888)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 10656 - EIF4A2
gene_id = 10656
de_isoform_probesets = c(2196658)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 4072 - ENO2
gene_id = 4072
de_isoform_probesets = c(2058838)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 1615 - EPB41L3
gene_id = 1615
de_isoform_probesets = c(277648,2347608,2347610,2347607,2347604,2347602,277653,277644,2347605,277639,2347609,2347606)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 10398 - FGD5
gene_id = 10398
de_isoform_probesets = c(2446882,1264346,2446884,2446885,1264286,2446883,2446890,1264300,1264390,2446888,2446875,2446876,1264370,1264271,
				 2446886,2446877,1264381,2446879,2446889,1264313,2446880,2446881,2446878,1264325,1264355,2446887,2446891,1264385,1264376)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 12008 - HHIP
gene_id = 12008
de_isoform_probesets = c(2463406,2463407,1461488,2463412,1461492,2463413,1461495,1461453,2463414,2463408,1461470,2463409,1461477,2463410,2463411,1461462,1461483,2463415)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 10935 - IL10RB
gene_id = 10935
de_isoform_probesets = c(1377842)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 739 - LAMA3
gene_id = 739
de_isoform_probesets = c(123032,123068,123103,123137,123170,123202,123233,123263,123292,123320,123347,123373,123398,123422,123445,
123467,123488,123508,123527,123545,123562,123578,123593,123607,123620,123632,123643,123653,123662,123670,123677,123683,123688,123692,
123695,123697,2334135,2334136,2334137,2334138,2334139,2334140,2334141,2334142,2334143,2334144,2334145,2334146,2334147,2334148,
2334149,2334150,2334151,2334152,2334153,2334154,2334155,2334156,2334157,2334158,2334159,2334160,2334161,2334162,2334163,2334164,
2334165,2334166,2334167,2334168,2334169,2334170,2334171,2334172)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 4647 - MLPH
gene_id = 4647
de_isoform_probesets = c(627821)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 2596 - PLCB4
gene_id = 2596
de_isoform_probesets = c(400639)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 7930 - PPP2R1B
gene_id = 7930
de_isoform_probesets = c(982489)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 3002 - RC74
gene_id = 3002
de_isoform_probesets = c(445760)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 16941 - RCC1
gene_id = 16941
de_isoform_probesets = c(1696636,1696624)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 9224 - SSBP2
gene_id = 9224
de_isoform_probesets = c(2167184,2167196)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 13665 - TPST1
gene_id = 13665
de_isoform_probesets = c(1579825)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])


######################################################################################
#ALEXA Gene: 4475 - UMPS
gene_id = 4475
de_isoform1_probesets = c(2381692,602470,602475)
de_isoform2_probesets = c(602471)
length(de_isoform1_probesets)
length(de_isoform2_probesets)
de_isoform1_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform1_probesets),]
de_isoform2_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform2_probesets),]

#Differential expression of isoform A
length(de_isoform1_data[,"mip_v_fur_DE_U"])
mean(de_isoform1_data[,"mip_v_fur_DE_U"])
sd(de_isoform1_data[,"mip_v_fur_DE_U"])

#Significant p-value?
mip101 = c(de_isoform1_data[,6], de_isoform1_data[,7], de_isoform1_data[,8])
mip5fu = c(de_isoform1_data[,9], de_isoform1_data[,10], de_isoform1_data[,11])
t.test(x=mip101, y=mip5fu, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

#Differential expression of isoform B
length(de_isoform2_data[,"mip_v_fur_DE_U"])
mean(de_isoform2_data[,"mip_v_fur_DE_U"])
sd(de_isoform2_data[,"mip_v_fur_DE_U"])

#Significant p-value?
mip101 = c(de_isoform2_data[,6], de_isoform2_data[,7], de_isoform2_data[,8])
mip5fu = c(de_isoform2_data[,9], de_isoform2_data[,10], de_isoform2_data[,11])
t.test(x=mip101, y=mip5fu, alternative = c("two.sided"), mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)




#mean and SD of expression of Isoform A in MIP101
mean(c(de_isoform1_data[,6], de_isoform1_data[,7], de_isoform1_data[,8]))
sd(c(de_isoform1_data[,6], de_isoform1_data[,7], de_isoform1_data[,8]))

#mean expression of Isoform B in MIP101
mean(c(de_isoform2_data[,6], de_isoform2_data[,7], de_isoform2_data[,8]))
sd(c(de_isoform2_data[,6], de_isoform2_data[,7], de_isoform2_data[,8]))


#mean and SD of expression of Isoform A in MIP/5FU
mean(c(de_isoform1_data[,9], de_isoform1_data[,10], de_isoform1_data[,11]))
sd(c(de_isoform1_data[,9], de_isoform1_data[,10], de_isoform1_data[,11]))

#mean expression of Isoform B in MIP/5FU
mean(c(de_isoform2_data[,9], de_isoform2_data[,10], de_isoform2_data[,11]))
sd(c(de_isoform2_data[,9], de_isoform2_data[,10], de_isoform2_data[,11]))



#mean ratio of Isoform A/B in MIP101
2^(mean(de_isoform1_data[,"mip_LOG2_I_U"]) - mean(mean(de_isoform2_data[,"mip_LOG2_I_U"])))

#sd of Isoform A/B in MIP101
sd(c(mean(de_isoform1_data[,6]) - mean(de_isoform2_data[,6]), mean(de_isoform1_data[,7]) - mean(de_isoform2_data[,7]), mean(de_isoform1_data[,8]) - mean(de_isoform2_data[,8])))

#mean ratio of Isoform A/B in MIP/5FU
2^(mean(de_isoform1_data[,"fur_LOG2_I_U"]) - mean(mean(de_isoform2_data[,"fur_LOG2_I_U"])))

#sd of Isoform A/B in MIP/5FU
sd(c(mean(de_isoform1_data[,9]) - mean(de_isoform2_data[,9]), mean(de_isoform1_data[,10]) - mean(de_isoform2_data[,10]), mean(de_isoform1_data[,11]) - mean(de_isoform2_data[,11])))



######################################################################################
#ALEXA Gene: 9470 - ZNF185
gene_id = 9470
de_isoform_probesets = c(1162625)
length(de_isoform_probesets)
de_isoform_data = exon_exp_data[which(exon_exp_data[,"ProbeSet_ID"] %in% de_isoform_probesets),]

de_gene_probesets = canonical_exp_data[which(canonical_exp_data[,"AlexaGene_ID"] == gene_id),]
de_other_data = de_gene_probesets[which(!(de_gene_probesets[,"ProbeSet_ID"] %in% de_isoform_probesets)),]

length(de_isoform_data[,"mip_v_fur_DE_U"])
mean(de_isoform_data[,"mip_v_fur_DE_U"])
sd(de_isoform_data[,"mip_v_fur_DE_U"])

length(de_other_data[,"mip_v_fur_DE_U"])
mean(de_other_data[,"mip_v_fur_DE_U"])
sd(de_other_data[,"mip_v_fur_DE_U"])














