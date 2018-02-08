#run edgeR robust
library(edgeR)
read.table(file="featureCounts", header=TRUE, row.names=("Geneid"))->x
x.genes<-x[,1:5] #featureCounts has gene data (start, stop, strand) as first 5 columns
x.counts<-x[,6:ncol(x)] #raw read counts 
group=c(1,1,1,1,2,2,2,1,1,2,2) #1’s are control group and 2's are test group

y<-DGEList(counts=x.counts,group=group,gene=x.genes) #x.genes is completely optional 
y <- DGEList(counts=x, group=group) #alt statement
y <- calcNormFactors(y) #normalize read counts
table(rowSums(cpm(y))>0)["TRUE"] #count number of expressed genes
table(rowSums(cpm(y))>1)["TRUE"] #count number of genes for some expression cutoff
keep<-rowSums(cpm(y))>1 #filter very low expressed genes
y<-y[keep,,keep.lib.sizes=FALSE] #filter step
y <- calcNormFactors(y) #re-normalize read counts with filtered genes

design <- model.matrix(~group)
y <- estimateGLMRobustDisp(y, design) #different robust step
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)

topTags(lrt, NROW(lrt$table)) -> geneList #outuput all genes
with(geneList$table, sum(FDR < 0.05)) #count DE genes
y.cpm<-cpm(y) #calculate cpm
y.rpkm<-rpkm(y) #calculate rpkm; if desired
final.output<-merge(geneList$table, y.cpm, by=0, all=TRUE) #merge tables by row names
write.table(final.output[order(final.output$FDR), ], file="edgeR-robust.tsv", sep="\t", row.names=F) #write results out
