#Written by Malachi Griffith
#Targeted examination of a single gene locus (INPP4B)

#Load neccessary libraries
library(gcrma)
library(preprocessCore)

#First load a junction expression matrix for known junctions only.
data_file1 = "/projects/alexa2/hmmSplicer/Summary/Compendium/JunctionExpressionMatrix.DA.tsv";
summary_file1 = "/projects/alexa2/hmmSplicer/Summary/Compendium/LibrarySummary.tsv";
data_file2 = "/projects/alexa2/hmmSplicer/Summary/SA_vs_Normals/JunctionExpressionMatrix.DA.tsv";
summary_file2 = "/projects/alexa2/hmmSplicer/Summary/SA_vs_Normals/LibrarySummary.tsv";
class_file2 = "/projects/alexa2/hmmSplicer/Summary/SA_vs_Normals/LibraryClasses.tsv"
data_file3 = "/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/JunctionExpressionMatrix.DA.tsv";
summary_file3 = "/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/LibrarySummary.tsv";
class_file3 = "/projects/alexa2/hmmSplicer/SA_TN_Breast/summary/LibraryClasses.tsv"

#Set output dir
setwd("/projects/alexa2/hmmSplicer/Summary/Compendium/INPP4B/")

#Set size of PDF output files
pdfh=12
pdfw=20


#Define those junctions that correspond to the long isoform only and those the could come from long or short (nothing is unique to the short isoform)
long=c("chr4:143545927-143569776(-)","chr4:143569820-143571772(-)","chr4:143571988-143603266(-)");
long_short=c("chr4:143169517-143222634(-)","chr4:143222788-143226747(-)","chr4:143226859-143248696(-)","chr4:143248793-143253145(-)","chr4:143253285-143262731(-)","chr4:143262848-143263895(-)","chr4:143264018-143265191(-)","chr4:143265363-143286443(-)","chr4:143286599-143300961(-)","chr4:143301164-143314235(-)","chr4:143314412-143333690(-)","chr4:143333798-143349028(-)","chr4:143349132-143349499(-)","chr4:143349629-143378467(-)","chr4:143378614-143401095(-)","chr4:143401167-143411266(-)","chr4:143411377-143446061(-)","chr4:143446140-143455315(-)","chr4:143455365-143543541(-)","chr4:143543657-143545809(-)");
total_j=c(long,long_short)

######################################################################################################################################################################
#Compendium data first
data = read.table(file=data_file1, header=TRUE, sep="\t", as.is=c(1,7,8,12), na.strings="na")
summary = read.table(file=summary_file1, header=TRUE, sep="\t", as.is=c(1), na.strings="na")

long_i=which(data[,"JID"]%in%long) 
long_short_i=which(data[,"JID"]%in%long_short)

#Create a normalized expression matrix to allow display of the levels from each library on the same scale - try quantiles normalized values
x = as.matrix(data[,13:(ncol(data)-2)])
x_norm = normalize.quantiles(x)

#Calculate expression levels for each isoform using the junctions defined above
#Long isoform will average of 3 junctions corresponding to first three exons
#Short isoform will be average of all other junctions MINUS the expression estimate of the long isoform
long_u = apply(x_norm[long_i,], 2, mean)
long_short_u = apply(x_norm[long_short_i,], 2, mean)

lib_ids = colnames(data[,13:(ncol(data)-2)])
inpp4b = data.frame(lib_ids, summary[,"Project"], long_u, long_short_u)
names(inpp4b) = c("Library_ID","Project","Long_U", "Long_Short_U")
inpp4b[,"Short_U"] = inpp4b[,"Long_Short_U"]-inpp4b[,"Long_U"]

#Reset negative values for the short isoform to 0 for display purposes...  these should be rare
inpp4b[which(inpp4b[,"Short_U"] < 0), "Short_U"] = 0

#Plot the expression of each isoform for each library - Use a line plot.  Define points according to project
pdf(file="INPP4B_Compendium.pdf", height=pdfh, width=pdfw)
par(font=2)
bp = barplot(height=t(as.matrix(inpp4b[,c("Short_U","Long_U")])), col=c("gray90","red"), border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2,
                      xlab="Library", ylab="Average Read Count per Junction per Library (for long/short isoforms)", main="Expression level of INPP4B isoforms (long+short) across 576 RNA-seq libraries",
                      names.arg=inpp4b[,"Project"], cex.names=0.6)
ltext=c("Short Isoform", "Long Isoform")
legend("topright", ltext, col=c("gray90","red"), pch=15)

cols=factor(inpp4b[,"Project"], label=1:12)
rcols=rainbow(12)
inpp4b[,"color_code"]=cols
inpp4b[,"color"]=rcols[inpp4b[,"color_code"]]
points(bp,rep(-1, length(bp)), pch=15, col=inpp4b[,"color"], cex=0.75)

ltext=levels(inpp4b[,"Project"])
legend("topleft", ltext, col=rainbow(12), pch=15)
dev.off()

######################################################################################################################################################################
#Now repeat the entire steps above but using a more limited set of Triple-Negative breast cancer samples plus other SA libraries and some normals
data = read.table(file=data_file2, header=TRUE, sep="\t", as.is=c(1,7,8,12), na.strings="na")
summary = read.table(file=summary_file2, header=TRUE, sep="\t", as.is=c(1), na.strings="na")
classes = read.table(file=class_file2, header=TRUE, sep="\t", as.is=c(1:8), na.strings="na")

#Define those junctions that correspond to the long isoform only and those the could come from long or short (nothing is unique to the short isoform)
long_i=which(data[,"JID"]%in%long) 
long_short_i=which(data[,"JID"]%in%long_short)

#Create a normalized expression matrix to allow display of the levels from each library on the same scale - try quantiles normalized values
x = as.matrix(data[,13:(ncol(data)-2)])
x_norm = normalize.quantiles(x)

#Calculate expression levels for each isoform using the junctions defined above
#Long isoform will average of 3 junctions corresponding to first three exons
#Short isoform will be average of all other junctions MINUS the expression estimate of the long isoform
long_u = apply(x_norm[long_i,], 2, mean)
long_short_u = apply(x_norm[long_short_i,], 2, mean)

lib_ids = colnames(data[,13:(ncol(data)-2)])
inpp4b = data.frame(lib_ids, summary[,"Project"], classes[,], long_u, long_short_u)
names(inpp4b) = c("Library_ID","Project", names(classes), "Long_U", "Long_Short_U")
inpp4b[,"Short_U"] = inpp4b[,"Long_Short_U"]-inpp4b[,"Long_U"]

#Reset negative values for the short isoform to 0 for display purposes...  these should be rare
inpp4b[which(inpp4b[,"Short_U"] < 0), "Short_U"] = 0

#Assign colors by 'project'
cols=factor(inpp4b[,"Project"], label=1:4)
rcols=rainbow(4)
inpp4b[,"color_code"]=cols
inpp4b[,"color"]=rcols[inpp4b[,"color_code"]]

#Flag 'normal' tissues
inpp4b[,"CancerNormal"] = "Cancer"
inpp4b[which(classes[,"Class_2"]=="Normal"),"CancerNormal"] = "Normal"

#Order according to Basal, vs. Non-Basal vs. TBD
inpp4b[,"BasalNonBasal"] = "TBD"
inpp4b[which(classes[,"Description"]=="Basal"),"BasalNonBasal"] = "Basal"
inpp4b[which(classes[,"Description"]=="Non-Basal"),"BasalNonBasal"] = "Non-Basal"
bas_o = order(inpp4b[,"BasalNonBasal"])
inpp4b_o = inpp4b[bas_o,]
normal_i = which(inpp4b_o[,"CancerNormal"]=="Normal")

basal_c = length(which(inpp4b_o[,"BasalNonBasal"]=="Basal"))
non_basal_c = length(which(inpp4b_o[,"BasalNonBasal"]=="Non-Basal"))
tbd_c = length(which(inpp4b_o[,"BasalNonBasal"]=="TBD"))
breast_c = length(which(inpp4b_o[,"Project"]=="SA_TN_Breast"))

#Plot the expression of each isoform for each library - Use a line plot.  Define points according to project
pdf(file="INPP4B_TN-Normals.pdf", height=pdfh, width=pdfw)
par(font=2)
bp = barplot(height=t(as.matrix(inpp4b_o[,c("Short_U","Long_U")])), col=c("gray90","red"), border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2,
                      xlab="Library", ylab="Average Read Count per Junction per Library (for long/short isoforms)", main="Expression level of INPP4B isoforms (long+short) across 152 RNA-seq libraries",
                      names.arg=inpp4b_o[,"Project"], cex.names=0.6)
ltext=c("Short Isoform", "Long Isoform")
legend("topright", ltext, col=c("gray90","red"), pch=15)

points(bp,rep(-1, length(bp)), pch=15, col=inpp4b_o[,"color"], cex=0.75)
points(bp[normal_i], inpp4b_o[normal_i,"Short_U"], pch=8, col="black", cex=2)

#Draw seperating lines for the basal, non-basal groups
offset = (bp[2]-bp[1])/2
basal_p = bp[basal_c]+offset
non_basal_p = bp[non_basal_c]+basal_p+offset
breast_p = bp[breast_c]+offset
abline(v=basal_p, lty=2, lwd=2, col="black")
abline(v=non_basal_p, lty=2, lwd=2, col="black")
abline(v=breast_p, lty=2, lwd=2, col="black")

ltext=c(levels(inpp4b_o[,"Project"]), "Normals")
legend("topleft", ltext, col=c(rainbow(4), "black"), pch=c(15,15,15,15,8))
text(x=15, y=100, "Basals", font=2, cex=1.5)
text(x=57, y=100, "Non-Basals", font=2, cex=1.5)
text(x=110, y=100, "Other", font=2, cex=1.5)
text(x=165, y=100, "Normal", font=2, cex=1.5)
dev.off()

#Write out data values to a tsv file
write.table(file="INPP4B.IsoformLevels.tsv", inpp4b, sep="\t", quote=FALSE)


######################################################################################################################################################################
#Finally, include only the Triple-Negatives and distinguish between Basal, Non-Basal and TBD
data = read.table(file=data_file3, header=TRUE, sep="\t", as.is=c(1,7,8,12), na.strings="na")
summary = read.table(file=summary_file3, header=TRUE, sep="\t", as.is=c(1), na.strings="na")
classes = read.table(file=class_file3, header=TRUE, sep="\t", as.is=c(1:8), na.strings="na")

#Define those junctions that correspond to the long isoform only and those the could come from long or short (nothing is unique to the short isoform)
long_i=which(data[,"JID"]%in%long) 
long_short_i=which(data[,"JID"]%in%long_short)
total_j_i=which(data[,"JID"]%in%total_j)

#Create a normalized expression matrix to allow display of the levels from each library on the same scale - try quantiles normalized values
x = as.matrix(data[,13:(ncol(data)-2)])
x_norm = normalize.quantiles(x)

#Calculate expression levels for each isoform using the junctions defined above
#Long isoform will average of 3 junctions corresponding to first three exons
#Short isoform will be average of all other junctions MINUS the expression estimate of the long isoform
long_u = apply(x_norm[long_i,], 2, mean)
long_short_u = apply(x_norm[long_short_i,], 2, mean)
total_u = apply(x_norm[total_j_i,], 2, mean)

lib_ids = colnames(data[,13:(ncol(data)-2)])
inpp4b = data.frame(lib_ids, summary[,"Project"], classes[,], long_u, long_short_u, total_u)
names(inpp4b) = c("Library_ID","Project", names(classes), "Long_U", "Long_Short_U","Total_U")
inpp4b[,"Short_U"] = inpp4b[,"Long_Short_U"]-inpp4b[,"Long_U"]
inpp4b[,"Prop_Long_U"] = (inpp4b[,"Long_U"]/inpp4b[,"Short_U"])


#Reset negative values for the short isoform to 0 for display purposes...  these should be rare
inpp4b[which(inpp4b[,"Short_U"] < 0), "Short_U"] = 0
inpp4b[which(inpp4b[,"Prop_Long_U"] < 0), "Prop_Long_U"] = 0

#Assign colors by 'Description'
cols=factor(inpp4b[,"Description"], label=1:2)
rcols=c("blue","green")
inpp4b[,"color_code"]=cols
inpp4b[,"color"]=rcols[inpp4b[,"color_code"]]

#Order according to Basal, vs. Non-Basal vs. TBD
bas_o = order(inpp4b[,"Description"])
inpp4b_o = inpp4b[bas_o,]

basal_c = length(which(inpp4b_o[,"Description"]=="Basal"))
non_basal_c = length(which(inpp4b_o[,"Description"]=="Non-Basal"))

#Plot the expression of each isoform for each library - Use a line plot.  Define points according to project
pdf(file="INPP4B_TripleNegative.pdf", height=pdfh, width=pdfw)
par(font=2)
bp = barplot(height=t(as.matrix(inpp4b_o[,c("Short_U","Long_U")])), col=c("gray90","red"), border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2,
                      xlab="Library", ylab="Average Read Count per Junction per Library (for long/short isoforms)", main="Expression level of INPP4B isoforms (long+short) across 84 RNA-seq libraries",
                      names.arg=inpp4b_o[,"Lbrary_Name"], cex.names=0.6, las=2)
ltext=c("Short Isoform", "Long Isoform")
legend("topright", ltext, col=c("gray90","red"), pch=15)

points(bp,rep(-1, length(bp)), pch=15, col=inpp4b_o[,"color"], cex=3)


#Draw seperating lines for the basal, non-basal groups
offset = (bp[2]-bp[1])/2
basal_p = bp[basal_c]+offset
non_basal_p = bp[non_basal_c]+basal_p+offset
abline(v=basal_p, lty=2, lwd=2, col="black")

ltext=names(table(inpp4b_o[,"Description"]))
legend("topleft", ltext, col=rcols, pch=15)
text(x=25, y=80, "Basals", font=2, cex=1.5)
text(x=75, y=80, "Non-Basals", font=2, cex=1.5)
dev.off()


#Try plotting the ratio of long/long+short
pdf(file="INPP4B_TripleNegative_Ratios.pdf", height=pdfh, width=pdfw)
par(font=2)
bp = barplot(height=log2(inpp4b_o[,"Prop_Long_U"]+1), col=c("light blue"), border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, ylim=c(-0.1,8),
                      xlab="Library", ylab="Proportion of long isoform expression (long/long+short)", main="Expression ratio of INPP4B isoforms (long/long+short) across 84 RNA-seq libraries",
                      names.arg=inpp4b_o[,"Lbrary_Name"], cex.names=0.6, las=2)
points(bp,rep(-0.1, length(bp)), pch=15, col=inpp4b_o[,"color"], cex=3)
offset = (bp[2]-bp[1])/2
basal_p = bp[basal_c]+offset
non_basal_p = bp[non_basal_c]+basal_p+offset
abline(v=basal_p, lty=2, lwd=2, col="black")
abline(v=non_basal_p, lty=2, lwd=2, col="black")

ltext=names(table(inpp4b_o[,"Description"]))
legend("topleft", ltext, col=rcols, pch=15)
text(x=25, y=80, "Basals", font=2, cex=1.5)
text(x=75, y=80, "Non-Basals", font=2, cex=1.5)
dev.off()


#Try some boxplots of the isoform levels for basal and non-basal libraries
pdf(file="INPP4B_TripleNegative_IsoformDist.pdf", height=7.5, width=7.5)
x=list(
inpp4b_o[which(inpp4b_o[,"Description"]=="Basal"),"Short_U"],
inpp4b_o[which(inpp4b_o[,"Description"]=="Non-Basal"),"Short_U"],
inpp4b_o[which(inpp4b_o[,"Description"]=="Basal"),"Long_U"],
inpp4b_o[which(inpp4b_o[,"Description"]=="Non-Basal"),"Long_U"]
)
par(font=2, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, font.sub=2)
names(x) = c("Basal", "Non-Basal", "Basal", "Non-Basal")
boxplot(x, col=c("gray90","gray90","red","red"), ylab="Average Read Count per Junction", xlab="Isoform & Type", main="INPP4B Isoform expression - Basal vs. Non-Basal TN Breast Cancers")
abline(v=2.5, lwd=2, lty=1)
text(x=1.2, y=60, "Short Isoform", font=2, cex=1.5)
text(x=3.5, y=60, "Long Isoform", font=2, cex=1.5)
dev.off()




