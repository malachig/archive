# load PAM50 data

pam50 <- read.delim("PAM50.centroids_unx.txt",sep="\t",header=T,row.names=6);
exprMat <-read.delim("TN_normalised_expr_gene_PAM50.txt",sep="\t",header=T,row.names=1);

M <- dim(exprMat)[[2]];
N <- dim(exprMat)[[1]];
IP <- order(row.names(pam50));
IE <- order(row.names(exprMat));

exprMat <- exprMat[IE,];
pam50 <- pam50[IP,];


exprClass <- matrix(nrow=5,ncol=M);
row.names(exprClass) <- names(pam50[1:5]);


# normalise each gene
scaled.mat <- as.matrix(t(scale(t(exprMat))));

# compute the spearman rank coefficient for each sample:class pair
for (m in 1:M) {
  for (k in 1:5) {
    exprClass[k,m] <- cor(scaled.mat[,m],pam50[,k],method="spearman");
  }
}

exprClass.df <- as.data.frame(exprClass);
names(exprClass.df) <- names(exprMat);

exprClass.df <- t(exprClass.df);

pam50Class <- vector(mode="character",length=M);
for (m in 1:M){
  pam50Class[m] <- names(pam50[1:5])[which.max(exprClass.df[m,])];
}

outDF <- cbind(exprClass.df,pam50Class);
names(outDF) <- c("Case","Basal","HER2","LuminalA","LuminalB","Normal","pam50Class");
write.table(outDF,file="PAM50_subtypes_TN.txt",sep="\t",row.names=T, col.names=T,quote=F);

                                                           
  



            
