#Define list of libraries of interest
libs_I=c("A00036", "A00042","A00053", "A00090", "A00091", "A00123", "A00126", "A00134", "A00155")

#Load AML library summary data file
j_file1="/projects/alexa2/hmmSplicer/Summary/AML_vs_Normals2/JunctionExpressionMatrix.tsv"
l_file1="/projects/alexa2/hmmSplicer/Summary/AML_vs_Normals2/LibrarySummary.tsv"
c_file1="/projects/alexa2/hmmSplicer/Summary/AML_vs_Normals2/LibraryClasses.tsv"

#Set the output path for figures
setwd("/projects/alexa2/hmmSplicer/Summary/AML_vs_Normals2/")
data=read.table(l_file1, header=TRUE, na.strings="na", as.is=c(1), sep="\t")
classes=read.table(c_file1, header=TRUE, na.strings="na", as.is=c(1:8), sep="\t")

#Set size of PDF output files
pdfh=12
pdfw=20


#####################################################################################################################################################################################
#First consider junctions by class (i.e. DA, NDA, D, A, N)
#For each category, calculate the proportion of the library that consists on that junction type versus all junctions.

#Determine the sum of each junction type
data[,"Total_1"]=apply(data[,c("DA_1","NDA_1","D_1","A_1","N_1")], 1, sum)
data[,"Total_1"]=apply(data[,c("DA_10","NDA_10","D_10","A_10","N_10")], 1, sum)



#Get the percent of each type (1 or more reads)
data[,"DAp"] =(data[,"DA_1"]/data[,"Total_1"])*100
data[,"NDAp"] =(data[,"NDA_1"]/data[,"Total_1"])*100
data[,"Dp"] =(data[,"D_1"]/data[,"Total_1"])*100
data[,"Ap"] =(data[,"A_1"]/data[,"Total_1"])*100
data[,"Np"] =(data[,"N_1"]/data[,"Total_1"])*100

data[,"DAp"] =(data[,"DA_10"]/data[,"Total_1"])*100
data[,"NDAp"] =(data[,"NDA_10"]/data[,"Total_1"])*100
data[,"Dp"] =(data[,"D_10"]/data[,"Total_1"])*100
data[,"Ap"] =(data[,"A_10"]/data[,"Total_1"])*100
data[,"Np"] =(data[,"N_10"]/data[,"Total_1"])*100

#Get the order according to each junction type
o_SIZE=order(data[,"Junction_Reads"])
o_DA=order(data[,"DAp"])
o_NDA=order(data[,"NDAp"])
o_D=order(data[,"Dp"])
o_A=order(data[,"Ap"])
o_N=order(data[,"Np"])

#Define colors according to the project
cols=as.character(factor(data[,"Project"], label=rainbow(4)))

#Libsize
pdf(file="AML_plus_Normals_LibSizes.pdf", height=pdfh, width=pdfw)
bar = barplot(height=data[o_SIZE,c("Junction_Reads")], col=cols[o_SIZE], border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2,
                      xlab="Library", ylab="Total junction counts", main="Library Sizes for 202 RNA-seq libraries",
                      names.arg=data[o_SIZE,"Project"], cex.names=0.6)

ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))

#Flag the splice factor mutated samples
z = which(data[o_SIZE,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_SIZE[z],c("Junction_Reads")], pch=8, col="black", cex=2)
dev.off()


#NDA
pdf("AML_NovelJunctionsByLibrary.pdf", height=pdfh, width=pdfw)
bar=barplot(height=data[o_NDA,c("NDAp")], col=cols[o_NDA], xlab="AML Library", ylab="% novel junctions ('NDA') (of total junctions expressed)", main="202 RNA-seq libraries (novel junctions)",
            border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, cex.names=0.6)
ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))
med_pos=length(o_NDA)/2
abline(v=bar[med_pos], lwd=2)

#Flag the splice factor mutated samples
z = which(data[o_NDA,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_NDA[z],c("NDAp")], pch=8, col="black", cex=2)
dev.off()


#####################################################################################################################################################################################
#Next consider the proportion of junctions that involve skipping of 1 or more exons
#This will be broken down into known and novel

#Determine the sum of each junction skipping type
#data[,"Total_2"]=apply(data[,c("ES_Known_0s", "ES_Known_1plus", "ES_Novel_0s", "ES_Novel_1plus")], 1, sum)
data[,"Total_2"]=data[,"Distinct_Junctions_0r"]



#Get the percent of each type of exon skipping
data[,"ES_Known_0s_p"] =(data[,"ES_Known_0s"]/data[,"Total_2"])*100
data[,"ES_Known_1s_p"] =(data[,"ES_Known_1s"]/data[,"Total_2"])*100
data[,"ES_Known_1plus_p"] =(data[,"ES_Known_1plus"]/data[,"Total_2"])*100
data[,"ES_Novel_0s_p"] =(data[,"ES_Novel_0s"]/data[,"Total_2"])*100
data[,"ES_Novel_1s_p"] =(data[,"ES_Novel_1s"]/data[,"Total_2"])*100
data[,"ES_Novel_1plus_p"] =(data[,"ES_Novel_1plus"]/data[,"Total_2"])*100

#Get the order according to each junction skipping type
o_Known_0s = order(data[,"ES_Known_0s_p"])
o_Known_1s = order(data[,"ES_Known_1s_p"])
o_Known_1plus = order(data[,"ES_Known_1plus_p"])
o_Novel_0s = order(data[,"ES_Novel_0s_p"])
o_Novel_1s = order(data[,"ES_Novel_1s_p"])
o_Novel_1plus = order(data[,"ES_Novel_1plus_p"])



#Known exon skipping (1 exon skipped only)
pdf("AML_KnownSkips_1s.pdf", height=pdfh, width=pdfw)
bar=barplot(height=data[o_Known_1s,c("ES_Known_1s_p")], col=cols[o_Known_1s], xlab="AML Library", ylab="% known single exon skipping (of total junctions expressed)", main="202 RNA-seq libraries (known single exon skips)",
            border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, cex.names=0.6)
ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))
med_pos=length(o_Known_1s)/2
abline(v=bar[med_pos], lwd=2)

#Flag the splice factor mutated samples
z = which(data[o_Known_1s,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_Known_1s[z],c("ES_Known_1s_p")], pch=8, col="black", cex=2)
dev.off()

#Known exon skipping (1+ exon skipped)
pdf("AML_KnownSkips_1plus.pdf", height=pdfh, width=pdfw)
bar=barplot(height=data[o_Known_1plus,c("ES_Known_1plus_p")], col=cols[o_Known_1plus], xlab="AML Library", ylab="% known exon skipping (of total junctions expressed)", main="202 RNA-seq libraries (known exon skips)",
            border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, cex.names=0.6)
ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))
med_pos=length(o_Known_1plus)/2
abline(v=bar[med_pos], lwd=2)

#Flag the splice factor mutated samples
z = which(data[o_Known_1plus,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_Known_1plus[z],c("ES_Known_1plus_p")], pch=8, col="black", cex=2)
dev.off()


#Novel exon skipping (1 exon skipped only)
pdf("AML_NovelSkips_1s.pdf", height=pdfh, width=pdfw)
bar=barplot(height=data[o_Novel_1s,c("ES_Novel_1s_p")], col=cols[o_Novel_1s], xlab="AML Library", ylab="% novel single exon skipping (of total junctions expressed)", main="202 RNA-seq libraries (novel single exon skips)",
            border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, cex.names=0.6)
ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))
med_pos=length(o_Novel_1s)/2
abline(v=bar[med_pos], lwd=2)

#Flag the splice factor mutated samples
z = which(data[o_Novel_1s,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_Novel_1s[z],c("ES_Novel_1s_p")], pch=8, col="black", cex=2)
dev.off()


#Novel exon skipping (1+ exon skipped)
pdf("AML_NovelSkips_1plus.pdf", height=pdfh, width=pdfw)
bar=barplot(height=data[o_Novel_1plus,c("ES_Novel_1plus_p")], col=cols[o_Novel_1plus], xlab="AML Library", ylab="% novel exon skipping (of total junctions expressed)", main="202 RNA-seq libraries (novel exon skips)",
            border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, cex.names=0.6)
ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))
med_pos=length(o_Novel_1plus)/2
abline(v=bar[med_pos], lwd=2)

#Flag the splice factor mutated samples
z = which(data[o_Novel_1plus,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_Novel_1plus[z],c("ES_Novel_1plus_p")], pch=8, col="black", cex=2)
dev.off()




#Get a library size normalized estimate of total exon skipping per library
data[,"ES_Any_Norm"] = ((data[,"ES_Known_1plus"]+data[,"ES_Novel_1plus"])/data[,"Junction_Reads"])*1000000
o_ES_Any_Norm = order(data[,"ES_Any_Norm"])

data[,"ES_Novel_Norm"] = ((data[,"ES_Novel_1plus"])/data[,"Junction_Reads"])*1000000
o_ES_Novel_Norm = order(data[,"ES_Novel_Norm"])

#Novel+Norm ('Any')
pdf("AML_AnySkips_1plus_Norm.pdf", height=pdfh, width=pdfw)
bar=barplot(height=data[o_ES_Any_Norm,c("ES_Any_Norm")], col=cols[o_ES_Any_Norm], xlab="AML Library", ylab="exon skipping (scaled by library size)", main="202 RNA-seq libraries (known+novel exon skips)",
            border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, cex.names=0.6, ylim=c(0,max(data[,c("ES_Any_Norm")])+100))
ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))
med_pos=length(o_ES_Any_Norm)/2
abline(v=bar[med_pos], lwd=2)

#Flag the splice factor mutated samples
z = which(data[o_ES_Any_Norm,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_ES_Any_Norm[z],c("ES_Any_Norm")], pch=8, col="black", cex=2)
dev.off();

#Novel only
pdf("AML_NovelSkips_1plus_Norm.pdf", height=pdfh, width=pdfw)
bar=barplot(height=data[o_ES_Novel_Norm,c("ES_Novel_Norm")], col=cols[o_ES_Novel_Norm], xlab="AML Library", ylab="exon skipping (scaled by library size)", main="202 RNA-seq libraries (novel exon skips)",
            border=NA, font.main = 2, font.lab = 2, font.axis=2, font.lab=2, cex.names=0.6, ylim=c(0,max(data[,c("ES_Novel_Norm")])+100))
ltext=c(levels(data[,"Project"]),"Mutated")
legend("topleft", ltext, col=c(rainbow(4),"black"), pch=c(rep(15,4),8))
med_pos=length(o_ES_Novel_Norm)/2
abline(v=bar[med_pos], lwd=2)

#Flag the splice factor mutated samples
z = which(data[o_ES_Novel_Norm,"Library_ID"] %in% libs_I)
points(x=bar[z], y=data[o_ES_Novel_Norm[z],c("ES_Novel_Norm")], pch=8, col="black", cex=2)
dev.off();

temp=data.frame(data[,c("Library_ID","Junction_Reads","Distinct_Junctions_0r","ES_Any_Norm","ES_Novel_Norm")], classes[,c("Lbrary_Name","Description")])
write.table(temp, file="temp.txt", sep="\t", quote=FALSE)

