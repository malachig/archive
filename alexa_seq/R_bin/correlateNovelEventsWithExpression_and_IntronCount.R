#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Are the number of novel splicing events observed for a gene correlated with the number of introns or expression level of the gene
#If these events are purely the function of noise from the splicing machinery then there should be a correlation as reported by Melamud and Moult (2009)

dir = "/projects/malachig/solexa/figures_and_stats/temp"
setwd(dir)

data = read.table("Gene_Expression_IntronCount_NovelJunction.txt", header=FALSE, sep=" ", as.is=1)
names(data) = c("ensg", "expression", "introns", "novel_junctions")

cor(x=data[,"expression"], y=data[,"novel_junctions"], method="spearman")
cor(x=data[,"introns"], y=data[,"novel_junctions"], method="spearman")


