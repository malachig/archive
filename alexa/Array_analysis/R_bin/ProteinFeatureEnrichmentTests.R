#Get data for simple Chi-Squared or Fisher's Exact test of Significant vs. Non-significant groups and categories (probeset position, etc.)

datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/"
dir(datadir)
setwd(datadir)

temp = read.table(file="Probeset_List_ANNOTATED.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(3,10,11,12), comment.char="")
complete_probeset_list = temp[which(temp[,"Probe_Type"]!="Control-Negative"),]

#SIG DE LIST
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted"
dir(datadir)
setwd(datadir)
sig_probeset_list = read.table(file="SigEvents_Log2Exp_FilteredOnExp_DE_ANNOTATED.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(3,13,14,15), comment.char="")

#SIG SI LIST
#datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/SI_values/MTP_formatted/"
#dir(datadir)
#setwd(datadir)
#sig_probeset_list = read.table(file="SigEvents_Log2Exp_FilteredOnExp_SI_ANNOTATED.txt", header=TRUE, sep="\t", na.strings = c("NA","na"), as.is=c(3,13,14,15), comment.char="")


#PEARSON CHI-SQUARED TEST and FISHER EXACT TESTS
#Generally it seems that you can use a Fisher Exact test for small samples and larger samples to a point.
#Eventually the sample size gets too large for computation of Fisher Exact, but with large counts the Chi-Squared test should be fine

#1.) Distribution of probe types
table(complete_probeset_list[,"Probe_Type"])
complete_E_count = length(which(complete_probeset_list[,"Probe_Type"]=="Exon"))
complete_C_count = length(which(complete_probeset_list[,"Probe_Type"]=="Exon-Exon" & complete_probeset_list[,"Exons_Skipped"] == 0))
complete_S_count = length(which(complete_probeset_list[,"Probe_Type"]=="Exon-Exon" & complete_probeset_list[,"Exons_Skipped"] > 0))
complete_EB_count = length(which(complete_probeset_list[,"Probe_Type"]=="Exon-Intron" | complete_probeset_list[,"Probe_Type"]=="Intron-Exon"))
complete_total = complete_E_count + complete_C_count + complete_S_count + complete_EB_count
(complete_E_count/complete_total)*100
(complete_C_count/complete_total)*100
(complete_S_count/complete_total)*100
(complete_EB_count/complete_total)*100

table(sig_probeset_list[,"Probe_Type"])
sig_E_count = length(which(sig_probeset_list[,"Probe_Type"]=="Exon"))
sig_C_count = length(which(sig_probeset_list[,"Probe_Type"]=="Exon-Exon" & sig_probeset_list[,"Exons_Skipped"] == 0))
sig_S_count = length(which(sig_probeset_list[,"Probe_Type"]=="Exon-Exon" & sig_probeset_list[,"Exons_Skipped"] > 0))
sig_EB_count = length(which(sig_probeset_list[,"Probe_Type"]=="Exon-Intron" | sig_probeset_list[,"Probe_Type"]=="Intron-Exon"))
sig_total = sig_E_count + sig_C_count + sig_S_count + sig_EB_count
(sig_E_count/sig_total)*100
(sig_C_count/sig_total)*100
(sig_S_count/sig_total)*100
(sig_EB_count/sig_total)*100

probe_types = matrix(c(complete_E_count, complete_C_count, complete_S_count, complete_EB_count, sig_E_count, sig_C_count, sig_S_count, sig_EB_count), nc = 2)
chisq.test(probe_types, y = NULL, correct = TRUE, p = rep(1/length(x), length(x)), rescale.p = FALSE, 
	     simulate.p.value = FALSE, B = 2000)


#2.) Distribution of probe positions
length(complete_probeset_list[,1])
length(sig_probeset_list[,1])
z1 = table(complete_probeset_list[,"ProbePosition"])
z2 = table(sig_probeset_list[,"ProbePosition"])
(z2["5prime_UTR"] / (z2["5prime_UTR"] + z2["ORF"] + z2["3prime_UTR"]))*100
(z2["ORF"] / (z2["5prime_UTR"] + z2["ORF"] + z2["3prime_UTR"]))*100
(z2["3prime_UTR"] / (z2["5prime_UTR"] + z2["ORF"] + z2["3prime_UTR"]))*100

probe_positions = matrix(c(z1["3prime_UTR"], z1["5prime_UTR"], z1["ORF"], z2["3prime_UTR"], z2["5prime_UTR"], z2["ORF"]), nc = 2)
chisq.test(probe_positions, y = NULL, correct = TRUE, p = rep(1/length(x), length(x)), rescale.p = FALSE, 
	     simulate.p.value = FALSE, B = 2000)

fisher.test(probe_positions, workspace = 2000000, hybrid = FALSE, control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)



#3.) Now consider the number of each type of protein feature in the group identified as significant compared to all probesets
(length(which(complete_probeset_list[,"TMD_count"] > 0)) / length(complete_probeset_list[,1]))*100
(length(which(complete_probeset_list[,"SignalPeptide_count"] > 0)) / length(complete_probeset_list[,1]))*100
(length(which(complete_probeset_list[,"CoiledCoil_count"] > 0)) / length(complete_probeset_list[,1]))*100
(length(which(complete_probeset_list[,"ProteinDomainCount"] > 0)) / length(complete_probeset_list[,1]))*100

(length(which(sig_probeset_list[,"TMD_count"] > 0)) / length(sig_probeset_list[,1]))*100
(length(which(sig_probeset_list[,"SignalPeptide_count"] > 0)) / length(sig_probeset_list[,1]))*100
(length(which(sig_probeset_list[,"CoiledCoil_count"] > 0)) / length(sig_probeset_list[,1]))*100
(length(which(sig_probeset_list[,"ProteinDomainCount"] > 0)) / length(sig_probeset_list[,1]))*100

#TMDs
tmd_complete = length(which(complete_probeset_list[,"TMD_count"] > 0))
no_tmd_complete = length(which(complete_probeset_list[,"TMD_count"] == 0))
tmd_sig = length(which(sig_probeset_list[,"TMD_count"] > 0))
no_tmd_sig = length(which(sig_probeset_list[,"TMD_count"] == 0))
tmds = matrix(c(tmd_complete, no_tmd_complete, tmd_sig, no_tmd_sig), nr = 2)
tmds
fisher.test(tmds, workspace = 200000, hybrid = FALSE, control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)

#SPs
sp_complete = length(which(complete_probeset_list[,"SignalPeptide_count"] > 0))
no_sp_complete = length(which(complete_probeset_list[,"SignalPeptide_count"] == 0))
sp_sig = length(which(sig_probeset_list[,"SignalPeptide_count"] > 0))
no_sp_sig = length(which(sig_probeset_list[,"SignalPeptide_count"] == 0))
sps <- matrix(c(sp_complete, no_sp_complete, sp_sig, no_sp_sig), nr = 2)
sps
fisher.test(sps, workspace = 200000, hybrid = FALSE, control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)

#CCs
cc_complete = length(which(complete_probeset_list[,"CoiledCoil_count"] > 0))
no_cc_complete = length(which(complete_probeset_list[,"CoiledCoil_count"] == 0))
cc_sig = length(which(sig_probeset_list[,"CoiledCoil_count"] > 0))
no_cc_sig = length(which(sig_probeset_list[,"CoiledCoil_count"] == 0))
ccs <- matrix(c(cc_complete, no_cc_complete, cc_sig, no_cc_sig), nr = 2)
ccs
fisher.test(ccs, workspace = 200000, hybrid = FALSE, control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)

#PDs
pd_complete = length(which(complete_probeset_list[,"ProteinDomainCount"] > 0))
no_pd_complete = length(which(complete_probeset_list[,"ProteinDomainCount"] == 0))
pd_sig = length(which(sig_probeset_list[,"ProteinDomainCount"] > 0))
no_pd_sig = length(which(sig_probeset_list[,"ProteinDomainCount"] == 0))
pds <- matrix(c(pd_complete, no_pd_complete, pd_sig, no_pd_sig), nr = 2)
pds
fisher.test(pds, workspace = 200000, hybrid = FALSE, control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)

#4.) For EST and mRNA support tests, only consider the 'alternative' probesets
EI_indexes = which(complete_probeset_list[,"Probe_Type"]=="Exon-Intron")
IE_indexes = which(complete_probeset_list[,"Probe_Type"]=="Intron-Exon")
SKIP_indexes = which(complete_probeset_list[,"Probe_Type"]=="Exon-Exon" & complete_probeset_list[,"Exons_Skipped"] > 0)
ALT_indexes = c(SKIP_indexes,EI_indexes,IE_indexes)
complete_alternative_probeset_list = complete_probeset_list[ALT_indexes,]

EI_indexes = which(sig_probeset_list[,"Probe_Type"]=="Exon-Intron")
IE_indexes = which(sig_probeset_list[,"Probe_Type"]=="Intron-Exon")
SKIP_indexes = which(sig_probeset_list[,"Probe_Type"]=="Exon-Exon" & sig_probeset_list[,"Exons_Skipped"] > 0)
ALT_indexes = c(SKIP_indexes,EI_indexes,IE_indexes)
sig_alternative_probeset_list = sig_probeset_list[ALT_indexes,]

#mRNA Evidence
known_complete = length(which(complete_alternative_probeset_list[,"Novel_or_Known"] == "known"))
novel_complete = length(which(complete_alternative_probeset_list[,"Novel_or_Known"] == "novel"))
known_sig = length(which(sig_alternative_probeset_list[,"Novel_or_Known"] == "known"))
novel_sig = length(which(sig_alternative_probeset_list[,"Novel_or_Known"] == "novel"))
mrna <- matrix(c(known_complete, novel_complete, known_sig, novel_sig), nr = 2)
mrna
fisher.test(mrna, workspace = 200000, hybrid = FALSE, control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)

#EST Support
est_complete = length(which(complete_alternative_probeset_list[,"EST_Support"] == "YES"))
no_est_complete = length(which(complete_alternative_probeset_list[,"EST_Support"] == "NO"))
est_sig = length(which(sig_alternative_probeset_list[,"EST_Support"] == "YES"))
no_est_sig = length(which(sig_alternative_probeset_list[,"EST_Support"] == "NO"))
est <- matrix(c(est_complete, no_est_complete, est_sig, no_est_sig), nr = 2)
est
fisher.test(est, workspace = 200000, hybrid = FALSE, control = list(), or = 1, alternative = "two.sided",
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, B = 2000)



#Permutation test
#Randomly sample populations of the same size as the 'significant' population with a similar composition of probe types
#Testing for enrichement of significant events in UTRs compared to ORF
#Determine how often you get as many 5' UTRs, as many 3' UTRs and as few ORFs
#Also determine how often you get as many probesets overlaping TMDs, SPs, CCs, ProteinDomains

sample_size = length(sig_probeset_list[,1])
probe_type_size = table(sig_probeset_list[,"Probe_Type"])

x = which(sig_probeset_list[,"Probe_Type"]=="Exon-Exon" & sig_probeset_list[,"Exons_Skipped"] > 0)
y1 = which(sig_probeset_list[,"Probe_Type"]=="Exon-Intron")
y2 = which(sig_probeset_list[,"Probe_Type"]=="Intron-Exon")
sig_alt_indexes = c(x,y1,y2)
probe_type_size_alt = length(sig_alt_indexes)

observed_positions = table(sig_probeset_list[,"ProbePosition"])
observed_known = length(which(sig_probeset_list[sig_alt_indexes,"Novel_or_Known"] == "known")) #only ALT probes considered
observed_est = length(which(sig_probeset_list[sig_alt_indexes,"EST_Support"] == "YES")) #only ALT probes considered
observed_tmd = length(which(sig_probeset_list[,"TMD_count"] > 0))
observed_sp = length(which(sig_probeset_list[,"SignalPeptide_count"] > 0))
observed_cc = length(which(sig_probeset_list[,"CoiledCoil_count"] > 0))
observed_pd = length(which(sig_probeset_list[,"ProteinDomainCount"] > 0))
 
iter = 1000
random_test = matrix(data = NA, nrow = iter, ncol = 9, dimnames = list(1:iter, c("five_prime", "orf", "three_prime", "tmd_count", "sp_count", "cc_count", "pd_count", "known_count", "est_count")))

E_indexes = which(complete_probeset_list[,"Probe_Type"]=="Exon")
EJ_indexes = which(complete_probeset_list[,"Probe_Type"]=="Exon-Exon")
EI_indexes = which(complete_probeset_list[,"Probe_Type"]=="Exon-Intron")
IE_indexes = which(complete_probeset_list[,"Probe_Type"]=="Intron-Exon")
I_indexes = which(complete_probeset_list[,"Probe_Type"]=="Intron")
SKIP_indexes = which(complete_probeset_list[,"Probe_Type"]=="Exon-Exon" & complete_probeset_list[,"Exons_Skipped"] > 0)
ALT_indexes = c(SKIP_indexes,EI_indexes,IE_indexes)

#Note that there were no significant Intron events, so drop this one
for (i in 1:iter){
  #Get a random sample of the p-values observed for each dataset
  rand_E_events = sample(E_indexes, probe_type_size["Exon"], replace = FALSE)
  rand_EJ_events = sample(EJ_indexes, probe_type_size["Exon-Exon"], replace = FALSE)
  rand_EI_events = sample(EI_indexes, probe_type_size["Exon-Intron"], replace = FALSE)
  rand_IE_events = sample(IE_indexes, probe_type_size["Intron-Exon"], replace = FALSE)
  rand_ALT_events = sample(ALT_indexes, probe_type_size_alt, replace = FALSE)

  #POSITIONS
  E_temp = complete_probeset_list[rand_E_events, ]
  EJ_temp = complete_probeset_list[rand_EJ_events, ]
  EI_temp = complete_probeset_list[rand_EI_events, ]
  IE_temp = complete_probeset_list[rand_IE_events, ]
  ALT_temp = complete_probeset_list[rand_ALT_events, ]

  E_pos_counts = table(E_temp[,"ProbePosition"])
  EJ_pos_counts = table(EJ_temp[,"ProbePosition"])
  EI_pos_counts = table(EI_temp[,"ProbePosition"])
  IE_pos_counts = table(IE_temp[,"ProbePosition"])

  five_prime_count = (E_pos_counts["5prime_UTR"] + EJ_pos_counts["5prime_UTR"] + EI_pos_counts["5prime_UTR"] + IE_pos_counts["5prime_UTR"])  
  orf_count = (E_pos_counts["ORF"] + EJ_pos_counts["ORF"] + EI_pos_counts["ORF"] + IE_pos_counts["ORF"])  
  three_prime_count = (E_pos_counts["3prime_UTR"] + EJ_pos_counts["3prime_UTR"] + EI_pos_counts["3prime_UTR"] + IE_pos_counts["3prime_UTR"])  
  
  random_test[i, "five_prime"] = five_prime_count
  random_test[i, "orf"] = orf_count
  random_test[i, "three_prime"] = three_prime_count

  #Known vs. Unknown? (Only consider skipping and boundary probesets for this test!!)
  ALT_known_count = length(which(ALT_temp[,"Novel_or_Known"] == "known"))
  random_test[i, "known_count"] = ALT_known_count
  cum_pvalue_known = (length(which(random_test[,"known_count"] >= observed_known))) / i

  #EST Support? (Only consider skipping and boundary probesets for this test!!)
  ALT_est_count = length(which(ALT_temp[,"EST_Support"] == "YES"))  
  random_test[i, "est_count"] = ALT_est_count
  cum_pvalue_est = (length(which(random_test[,"est_count"] >= observed_est))) / i

  #DOMAINS
  E_tmd_count = length(which(E_temp[,"TMD_count"] > 0))  
  EJ_tmd_count = length(which(EJ_temp[,"TMD_count"] > 0))  
  EI_tmd_count = length(which(EI_temp[,"TMD_count"] > 0))  
  IE_tmd_count = length(which(IE_temp[,"TMD_count"] > 0))  
  tmd_count = E_tmd_count + EJ_tmd_count + EI_tmd_count + IE_tmd_count
  random_test[i, "tmd_count"] = tmd_count

  E_sp_count = length(which(E_temp[,"SignalPeptide_count"] > 0))  
  EJ_sp_count = length(which(EJ_temp[,"SignalPeptide_count"] > 0))  
  EI_sp_count = length(which(EI_temp[,"SignalPeptide_count"] > 0))  
  IE_sp_count = length(which(IE_temp[,"SignalPeptide_count"] > 0))  
  sp_count = E_sp_count + EJ_sp_count + EI_sp_count + IE_sp_count
  random_test[i, "sp_count"] = sp_count

  E_cc_count = length(which(E_temp[,"CoiledCoil_count"] > 0))  
  EJ_cc_count = length(which(EJ_temp[,"CoiledCoil_count"] > 0))  
  EI_cc_count = length(which(EI_temp[,"CoiledCoil_count"] > 0))  
  IE_cc_count = length(which(IE_temp[,"CoiledCoil_count"] > 0))  
  cc_count = E_cc_count + EJ_cc_count + EI_cc_count + IE_cc_count
  random_test[i, "cc_count"] = cc_count

  E_pd_count = length(which(E_temp[,"ProteinDomainCount"] > 0))  
  EJ_pd_count = length(which(EJ_temp[,"ProteinDomainCount"] > 0))  
  EI_pd_count = length(which(EI_temp[,"ProteinDomainCount"] > 0))  
  IE_pd_count = length(which(IE_temp[,"ProteinDomainCount"] > 0))  
  pd_count = E_pd_count + EJ_pd_count + EI_pd_count + IE_pd_count
  random_test[i, "pd_count"] = pd_count

  message1 = paste(i, ": Known_OR = ", cum_pvalue_known, ", EST_OR = ", cum_pvalue_est)
  print (message1)
}

#Dump the resulting random observations to a file for later reference
write.table(random_test, file = "SigEvents_RandomSampling.txt", append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE, qmethod = c("escape", "double"))

(median(random_test[,"five_prime"]) / length(sig_probeset_list[,1]))*100
(median(random_test[,"orf"]) / length(sig_probeset_list[,1]))*100
(median(random_test[,"three_prime"]) / length(sig_probeset_list[,1]))*100
(median(random_test[,"known_count"]) / probe_type_size_alt)*100
(median(random_test[,"est_count"]) / probe_type_size_alt)*100
(median(random_test[,"tmd_count"]) / length(sig_probeset_list[,1]))*100
(median(random_test[,"sp_count"]) / length(sig_probeset_list[,1]))*100
(median(random_test[,"cc_count"]) / length(sig_probeset_list[,1]))*100
(median(random_test[,"pd_count"]) / length(sig_probeset_list[,1]))*100

#Final P-value calculations

#Enriched compared to random
(length(which(random_test[,"five_prime"] >= observed_positions["5prime_UTR"]))) / iter 
(length(which(random_test[,"orf"] >= observed_positions["ORF"]))) / iter 
(length(which(random_test[,"three_prime"] >= observed_positions["3prime_UTR"]))) / iter 
(length(which(random_test[,"known_count"] >= observed_known))) / iter
(length(which(random_test[,"est_count"] >= observed_est))) / iter
(length(which(random_test[,"tmd_count"] >= observed_tmd))) / iter
(length(which(random_test[,"sp_count"] >= observed_sp))) / iter
(length(which(random_test[,"cc_count"] >= observed_cc))) / iter
(length(which(random_test[,"pd_count"] >= observed_pd))) / iter

#Depleted compared to random
(length(which(random_test[,"five_prime"] <= observed_positions["5prime_UTR"]))) / iter 
(length(which(random_test[,"orf"] <= observed_positions["ORF"]))) / iter 
(length(which(random_test[,"three_prime"] <= observed_positions["3prime_UTR"]))) / iter 
(length(which(random_test[,"known_count"] <= observed_known))) / iter
(length(which(random_test[,"est_count"] <= observed_est))) / iter
(length(which(random_test[,"tmd_count"] <= observed_tmd))) / iter
(length(which(random_test[,"sp_count"] <= observed_sp))) / iter
(length(which(random_test[,"cc_count"] <= observed_cc))) / iter
(length(which(random_test[,"pd_count"] <= observed_pd))) / iter

















