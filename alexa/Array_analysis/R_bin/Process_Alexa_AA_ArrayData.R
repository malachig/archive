#Written by Malachi Griffith
#Pre-Process data for ALEXA-AFFY common probesets to correct for background specific to each experiment and normalize across arrays

#Steps for each platform
#- Import array data for negative control probes
#- Import array data for Alexa/Affy common probesets
#- Background correct these values by fitting an appropriate linear model to the negative control data (scale to prevent negative values)
#- Use quanitiles normalization
#- Print out background corrected, normalized data values


library(gcrma)

#A.) PROCESS ALEXA DATA

#Get the raw data for negative control probes
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Raw"
dir(datadir)
setwd(datadir)
data_complete = read.table("All_hybes_withProbeInfo.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(4,5,12,15), na.strings='na')
data_control_alexa = data_complete[which(data_complete[,"Probe_Type"]=="Control-Negative"),]
data_control_alexa = data_control_alexa[,c("Probe_ID","ProbeSet_ID","MIP_EF_A","MIP_EF_B","MIP_GH_A","FUR_EF_A","FUR_EF_B","FUR_GH_A","Probe_Tm")]

#Get the raw data for Alexa/Affy probesets
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Raw/AA_probesets"
dir(datadir)
setwd(datadir)
data_AA_alexa = read.table("AA_alexa_data.txt", header=T, quote="", sep="\t", comment.char="", as.is=c(5,6,15,18), na.strings='na')
data_AA_alexa = data_AA_alexa[,c("AA_probeset_id","AlexaGene_ID","Probe_ID","ProbeSet_ID","MIP_EF_A","MIP_EF_B","MIP_GH_A","FUR_EF_A","FUR_EF_B","FUR_GH_A","Probe_Tm","Probe_Type","Exons_Skipped")]

#Fit a loess model to the negative control data
B_MIP_EF_A = fitLoessModelToTmMediansBins(data_control_alexa, "MIP_EF_A", min(data_control_alexa[,"Probe_Tm"]), max(data_control_alexa[,"Probe_Tm"]), 0.1)
B_FUR_EF_A = fitLoessModelToTmMediansBins(data_control_alexa, "FUR_EF_A", min(data_control_alexa[,"Probe_Tm"]), max(data_control_alexa[,"Probe_Tm"]), 0.1)
B_MIP_EF_B = fitLoessModelToTmMediansBins(data_control_alexa, "MIP_EF_B", min(data_control_alexa[,"Probe_Tm"]), max(data_control_alexa[,"Probe_Tm"]), 0.1)
B_FUR_EF_B = fitLoessModelToTmMediansBins(data_control_alexa, "FUR_EF_B", min(data_control_alexa[,"Probe_Tm"]), max(data_control_alexa[,"Probe_Tm"]), 0.1)
B_MIP_GH_A = fitLoessModelToTmMediansBins(data_control_alexa, "MIP_GH_A", min(data_control_alexa[,"Probe_Tm"]), max(data_control_alexa[,"Probe_Tm"]), 0.1)
B_FUR_GH_A = fitLoessModelToTmMediansBins(data_control_alexa, "FUR_GH_A", min(data_control_alexa[,"Probe_Tm"]), max(data_control_alexa[,"Probe_Tm"]), 0.1)

#Calculate the background subtracted values for each dataset
memory.limit(size = 2000)
bgc_MIP_EF_A = generate_BGC_values(data_AA_alexa, "MIP_EF_A", B_MIP_EF_A)
bgc_FUR_EF_A = generate_BGC_values(data_AA_alexa, "FUR_EF_A", B_FUR_EF_A)
bgc_MIP_EF_B = generate_BGC_values(data_AA_alexa, "MIP_EF_B", B_MIP_EF_B)
bgc_FUR_EF_B = generate_BGC_values(data_AA_alexa, "FUR_EF_B", B_FUR_EF_B)
bgc_MIP_GH_A = generate_BGC_values(data_AA_alexa, "MIP_GH_A", B_MIP_GH_A)
bgc_FUR_GH_A = generate_BGC_values(data_AA_alexa, "FUR_GH_A", B_FUR_GH_A)

#Scale values to prevent any -ve values in any dataset
data_AA_alexa[,"MIP_EF_A_bgc"] = bgc_MIP_EF_A + 16
data_AA_alexa[,"MIP_EF_B_bgc"] = bgc_MIP_EF_B + 16
data_AA_alexa[,"MIP_GH_A_bgc"] = bgc_MIP_GH_A + 16
data_AA_alexa[,"FUR_EF_A_bgc"] = bgc_FUR_EF_A + 16
data_AA_alexa[,"FUR_EF_B_bgc"] = bgc_FUR_EF_B + 16
data_AA_alexa[,"FUR_GH_A_bgc"] = bgc_FUR_GH_A + 16

#Does this BGC procedure actually reduce variability across replicates?  Not really ??!
median(apply(data_AA_alexa[,c("MIP_EF_A","MIP_EF_B","MIP_GH_A")], 1, sd))
median(apply(data_AA_alexa[,c("MIP_EF_A_bgc","MIP_EF_B_bgc","MIP_GH_A_bgc")], 1, sd))

median(apply(data_AA_alexa[,c("FUR_EF_A","FUR_EF_B","FUR_GH_A")], 1, sd))
median(apply(data_AA_alexa[,c("FUR_EF_A_bgc","FUR_EF_B_bgc","FUR_GH_A_bgc")], 1, sd))

#Now Normalize the data
data_simple = data.frame(data_AA_alexa[,"MIP_EF_A_bgc"], data_AA_alexa[,"MIP_EF_B_bgc"], data_AA_alexa[,"MIP_GH_A_bgc"], 
				 data_AA_alexa[,"FUR_EF_A_bgc"], data_AA_alexa[,"FUR_EF_B_bgc"], data_AA_alexa[,"FUR_GH_A_bgc"])

data_names = c("MIP_EF_A", "MIP_EF_B", "MIP_GH_A", "FUR_EF_A", "FUR_EF_B", "FUR_GH_A")
names(data_simple) = data_names
x = as.matrix(data_simple)

data_simple_norm = normalize.quantiles(x)
data_simple_norm = as.data.frame(data_simple_norm)
dimnames(data_simple_norm)[[2]] = data_names


#Now build a new dataframe like the orginal one with BGC and Normalized values instead of raw values
new_data = data_AA_alexa[,c("AA_probeset_id","AlexaGene_ID","Probe_ID","ProbeSet_ID")]
new_data[,"MIP_EF_A"] = round(data_simple_norm[,"MIP_EF_A"], digits = 3)
new_data[,"MIP_EF_B"] = round(data_simple_norm[,"MIP_EF_B"], digits = 3)
new_data[,"MIP_GH_A"] = round(data_simple_norm[,"MIP_GH_A"], digits = 3)
new_data[,"FUR_EF_A"] = round(data_simple_norm[,"FUR_EF_A"], digits = 3)
new_data[,"FUR_EF_B"] = round(data_simple_norm[,"FUR_EF_B"], digits = 3)
new_data[,"FUR_GH_A"] = round(data_simple_norm[,"FUR_GH_A"], digits = 3)
new_data[,"Probe_Tm"] = data_AA_alexa[,"Probe_Tm"]
new_data[,"Probe_Type"] = data_AA_alexa[,"Probe_Type"]
new_data[,"Exons_Skipped"] = "na"

#Change to a new output directory and write out the BGC, normalized data file
datadir = "E:/My Documents/Grad Studies/Project Documents/ArrayAnalysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Quantiles_TmBGC/AA_probesets"
dir(datadir)
setwd(datadir)

write.table(new_data[,], file="AA_alexa_data_BGC_Norm.txt", eol="\n", sep="\t", row.names=FALSE, quote=FALSE, na = "na")



###################################################################################################
#Create a function to do this for any set of negative control probes with a range of Tm
#Provide an array containing Tm in the first column, followed by labeled columns of experimental I values
#Specify the experiment name to consider, specify the min_tm and max_tm possible
#Specify the bin_size to use
#Return an array of median intensity values for each bin along with the max and min Tm of each bin
##################################################################################################
fitLoessModelToTmMediansBins = function(tm_I_array,exp_name,min_tm, max_tm, bin_size){
  tm_diff=(abs(max_tm-min_tm))
  bin_count = round((tm_diff/bin_size), digits=0)

  Tm_start = min_tm
  nc_median_tm = matrix(0,nrow=bin_count,ncol=6, dimnames=list(1:bin_count, c("Tm_start","Tm_end","Tm_bin_centre","median","mean","bin_size")))

  for(i in 1:bin_count){
    x=which((tm_I_array[,"Probe_Tm"] > Tm_start) & (tm_I_array[,"Probe_Tm"] < (Tm_start+bin_size)))
    y = tm_I_array[x,exp_name]

    #Identify all the I values in this Tm bin that are NOT outliers
    #Where outliers are defined as values that are less than Q1-(IQR*.5) or more than Q3+(IQR*.5)
    q1 = quantile(y, 0.25)
    q3 = quantile(y, 0.75)
    iqr = (IQR(y))*0.5
    lower_limit = q1 - iqr
    upper_limit = q3 + iqr

    non_outliers = which((y > lower_limit) & (y < upper_limit))

    bin_I = y[non_outliers]

    med_I = median(bin_I)
    mean_I = mean(bin_I)

    Tm_bin_centre = Tm_start + ((abs((Tm_start+bin_size)-Tm_start))/2)

    nc_median_tm[i,"median"] = med_I
    nc_median_tm[i,"mean"] = mean_I
    nc_median_tm[i,"Tm_start"] = Tm_start
    nc_median_tm[i,"Tm_end"] = Tm_start+bin_size
    nc_median_tm[i,"Tm_bin_centre"] = Tm_bin_centre
    nc_median_tm[i,"bin_size"] = length(non_outliers)
    Tm_start = Tm_start+bin_size
  }

  title = paste(exp_name, "- NC Probes")

  #Fit a linear model to the Tm/probe intensity relationship
  x = nc_median_tm[,"Tm_bin_centre"]
  y = nc_median_tm[,"median"]

  #Linear Fit
  #lmfit = lm(y~x)
  #fitted.y = b[1] + (b[2]*x)

  #Quadratic Fit
  #lmfit = lm(y~x + I((x)^2))
  #b = lmfit$coef 
  #fitted.y = (b[1]) + (b[2] * x) + (b[3] * (x^2))

  #Exponential Fit
  #b = expfit(x, y, 1000)
  #fitted.y = b[1]*(1-(exp(b[2]*x)))

  #Exponential Fit 2
  #cc = coef(lm(log(y) ~ x))
  #b = nls(y ~ exp(a + b*x), start = list(a = cc[1], b = cc[2]))
  #fitted.y = fitted(nls(y ~ exp(a + b*x), start = list(a = cc[1], b = cc[2])))
  #fitted.y = (b[[1]])*(1-(exp((b[[2]])*x)))

  #Loess
  loess_fitted = predict(loess(y~x))
  loess_fit = loess(y~x, control=loess.control(surface="direct"))

  plot(x, y, col="blue", ylab="Raw Intensity", xlab="Probe Tm (C)", main=title, pch=16)
  lines(x, loess_fitted, col="red", lwd=2)

  return(loess_fit)
}
####################################################################################################


####################################################################################################
#Using a loess model fit, generate background subtracted values                                   #
####################################################################################################
generate_BGC_values = function(data, exp, model_fit){
  bgc_column_name = paste(exp, "_bgc", sep="")

  total_values = length(data[,"Probe_Tm"])
  message1 = paste("Processing:",total_values,"total probe intensity values")
  print(message1)

  #Use the model provide to calculate background estimates
  #bg_values = model_fit[1] + (model_fit[2]* data[,"Probe_Tm"])
  #bg_values = (model_fit[1]) + (model_fit[2] * data[,"Probe_Tm"]) + (model_fit[3] * ((data[,"Probe_Tm"])^2))

  bg_values = (predict(model_fit, data[,"Probe_Tm"], se=TRUE))$fit

  #Subtract the estimated background values
  subtracted_values = data[,exp] - bg_values
  min_after_subtracting = min(subtracted_values)

  #Note the minimum observed value after this background correction
  message2 = paste("Found a minimum value of",min_after_subtracting,"after background correcting")
  print(message2)

  #Count the number of -ve values that result from the background subtraction
  neg_val_count = length(which (subtracted_values < 0))
  message3 = paste("Found a total of",neg_val_count,"negative values after background correcting")
  print(message3)

  subtracted_values[which(subtracted_values < 0)] = 0

  return(subtracted_values)
}
