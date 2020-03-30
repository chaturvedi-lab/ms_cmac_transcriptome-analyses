library(edgeR)
library(limma)
library(Glimma)
library(gplots)
#library(org.Mm.eg.db)
library(RColorBrewer)

#read in the file with all the counts
seqdata<-read.table("allsamples_featurecounts_annot.txt", header=T)
seqdata<-matrix(scan("allsamples_featurecounts_annot.txt",n=49512*30,sep="\t"),nrow=49512,ncol=30,byrow=TRUE)
head(seqdata)
#read in the sample ids to get treatments and lines
samples<-read.table("sampleids.txt", header=T)
group<-interaction(samples$line, samples$host)

#format the data
#remove first four columns from seqdata to make count matrix
countdata<-seqdata[,-c(1:4,30)]
#look at the output
head(countdata) 

#store gene IDs as rownames
rownames(countdata)<-seqdata[,1]
#look at the output
head(countdata) 

#rename samples to have line and host in the names
newids<-interaction(samples$sample, samples$line, samples$host)
colnames(countdata)<-newids
head(countdata)
##filtering to remove lowly expressed genes using edgeR
#obtain CPMs
myCPM<-cpm(countdata)
#have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are XX genes that have TRUEs in all 25 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep) #4411 25

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
pdf("cpm.pdf", width=11, height=10)
par(mar=c(5,5,3,1))
plot(myCPM[,1],countdata[,1], cex.lab=1.5)
dev.off()

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
pdf("cpm_lowcounts.pdf", width=11, height=10)
par(mar=c(5,5,3,1))
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)
dev.off()

#convert counts to DGEList object
y <- DGEList(counts.keep)
# have a look at y
y

# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples

##Quality control
##Library sizes and distribution plot
y$samples$lib.size
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
pdf("librarysize_barplot.pdf", width=11, height=10)
par(mar=c(8,5,3,1))
barplot(y$samples$lib.size,names=colnames(y),las=2, cex.lab=1.3)
# Add a title to the plot
title("Barplot of library sizes", cex=2)
dev.off()

#Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. Next we’ll use box plots to check the distribution of the read counts on the log2 scale. We can use the cpm function to get log2 counts per million, which are corrected for the different library sizes. The cpm function also adds a small offset to avoid taking log of zero.
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
pdf("countdistr_boxplot.pdf", width=11, height=10)
par(mar=c(8,5,3,1))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,cex.lab=1.3)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue", lwd=2)
title("Boxplots of logCPMs (unnormalised)", cex=2)
dev.off()

##input for the barplots below
dat<-read.csv("deg_numbers.txt", header = F)
dat_t<-t(dat)
colnames(dat_t)<-dat_t[1,]
#reorder columns
dat_t<-dat_t[-1,]
colnames(dat_t)[4]<-"M1.M - L1.L"
dat_o<-subset(dat_t, select=c(7,9,4,8,6,2,3,5,1))
7,9,6,4,8,
library(RColorBrewer)
coul <- brewer.pal(5, "Set2")[c(3,4,5)]
pdf("barplot_degs.pdf", width = 15, height = 12)
par(mar=c(10,8,3,4) + 0.1)
barplot(dat_t[-1,], col=coul , border="white", xlab="", ylab="Number of genes", las=2, cex.lab=2, cex.names = 1.5)
legend("topleft", inset=.02, c("Upregulated","Downregulated"), fill=coul[c(2,3)],cex=2)
dev.off()

#Multidimensional scaling (MDS) plot
library(RColorBrewer)
colors<-brewer.pal(n=6, name="Dark2") #color by group
pch<-c(17,17,19,17,19,17,17,19,17,17,19,17,17,19,17,19,17,19,17,17,19,17,19,17,19) #17=M 19=L
#plot
pdf("mds_barplot.pdf", width = 15, height = 16)
par(mfrow=c(2,1))
par(xpd = T, pty='s')
plotMDS(y, col=colors[as.numeric(group)], pch=pch, cex=2,ylim=c(-2,2.5), cex.lab=2,cex.main=2, main = "(A) Multidimensional scaling plot ")
legend("topright", inset = c(- 0.42, 0),legend=levels(group)[-3], pch=c(19,19,17,17,17), col=colors[-3], ncol=2, cex=1.5)
#title("(A) Multidimensional scaling plot", adj = 0.5, cex=2)

coul <- brewer.pal(5, "Set2")[c(4,5)]
#put = 'm'
par(pty='m')
par(mar=c(7,8,4,4))
barplot(dat_o[-1,], col=coul , border="white", xlab="", ylab="Number of genes", las=2, cex.lab=2,cex.main=2, cex.names = 1,main="(B) Barplot for differentially expressed genes")
legend("topright", inset=.02, c("Upregulated","Downregulated"), fill=coul,cex=1.5)
#title("(B) Barplot for differentially expressed genes", adj = 0.5, cex=2)
dev.off()

####heat map
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
select_var <- names(sort(var_genes, decreasing=TRUE))[1:50]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm) #50 25
head(highly_variable_lcpm)
colnames(highly_variable_lcpm)
#get gene ids for top 50 degs
top50degid<-seqdata[which(rownames(highly_variable_lcpm) %in% seqdata$Geneid),][,30]

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[samples$host]
temp<-highly_variable_lcpm
rownames(temp)<-top50degid

pdf("heatmap_degs.pdf", width=10,height=10)
#par(mar=c(8,5,5,7))
heatmap.2(temp,margins= c(8,8),cexRow=0.6,cexCol=1,labCol=top50degid,col=rev(morecols(50)),trace="none",main="Top 50 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()


##differential expression with limma-voom
#Voom transformation and calculation of variance weights
mm <- model.matrix(~0 + group)
par(mfrow=c(1,1))
v <- voom(y, mm, plot = T)
dev.copy(pdf, "voom_trendplot.pdf")
dev.off()

#details of the plot above
#What is voom doing?
#1. Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the normalization factors we calculated earlier
#2. A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
#3. A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
#4. The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.

#We can repeat the box plots for the normalised data to compare to before normalisation. The expression values in v$E are already log2 values so we don’t need to log-transform.
pdf("normdat_compare_boxplot.pdf", width=15, height=10)
par(mfrow=c(1,2))
par(mar=c(8,5,3,1))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM", cex.main=2)
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue", lwd=2)
par(mar=c(8,5,3,1))
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM",cex.main=2)
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue", lwd=2)
dev.off()


# Testing for differential expression in limma
#Fit linear model using weighted least squares for each gene using lmFit
fit <- lmFit(v, mm)
head(coef(fit))

#compare between groups as "contrasts" of these fitted linear models
#specify which groups to compare (trying just lentil vs mung)
#lentil 
contr_l1lvl1rl <- makeContrasts(groupL1.L - groupL1R.L, levels = colnames(coef(fit)))
contr_l1lvl1m <- makeContrasts(groupL1.L - groupL1.M, levels = colnames(coef(fit)))
contr_l1lvl1rm <- makeContrasts(groupL1.L - groupL1R.M, levels = colnames(coef(fit)))
contr_l1lvm1m <- makeContrasts(groupL1.L - groupM1.M, levels = colnames(coef(fit)))
#mung
contr_m1mvl1m <- makeContrasts(groupM1.M - groupL1.M, levels = colnames(coef(fit)))
contr_m1mvl1rl <- makeContrasts(groupM1.M - groupL1R.L, levels = colnames(coef(fit)))
contr_m1mvl1rm <- makeContrasts(groupM1.M - groupL1R.M, levels = colnames(coef(fit)))


#estimate contrast for each gene
#lentil
fit_l1lvl1rl<- contrasts.fit(fit, contr_l1lvl1rl)
fit_l1lvl1m<- contrasts.fit(fit, contr_l1lvl1m)
fit_l1lvl1rm<- contrasts.fit(fit, contr_l1lvl1rm)
fit_l1lvm1m<- contrasts.fit(fit, contr_l1lvm1m)
#mung
fit_m1mvl1m<- contrasts.fit(fit, contr_m1mvl1m)
fit_m1mvl1rl<- contrasts.fit(fit, contr_m1mvl1rl)
fit_m1mvl1rm<- contrasts.fit(fit, contr_m1mvl1rm)

#Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error)
#lentil
fit_l1lvl1rl<- eBayes(fit_l1lvl1rl)
fit_l1lvl1m<- eBayes(fit_l1lvl1m)
fit_l1lvl1rm<- eBayes(fit_l1lvl1rm)
fit_l1lvm1m<- eBayes(fit_l1lvm1m)
#mung
fit_m1mvl1m<- eBayes(fit_m1mvl1m)
fit_m1mvl1rl<- eBayes(fit_m1mvl1rl)
fit_m1mvl1rm<- eBayes(fit_m1mvl1rm)

#What genes are most differentially expressed?
#top.table <- topTable(tmp, sort.by = "P", n = Inf)
#head(top.table, 20)
#lentil
toptab_l1lvl1rl<- topTable(fit_l1lvl1rl, sort.by = "P", n = Inf)
toptab_l1lvl1m<- topTable(fit_l1lvl1m, sort.by = "P", n = Inf)
toptab_l1lvl1rm<- topTable(fit_l1lvl1rm, sort.by = "P", n = Inf)
toptab_l1lvm1m<- topTable(fit_l1lvm1m, sort.by = "P", n = Inf)
#mung
toptab_m1mvl1m<- topTable(fit_m1mvl1m, sort.by = "P", n = Inf)
toptab_m1mvl1rl<- topTable(fit_m1mvl1rl, sort.by = "P", n = Inf)
toptab_m1mvl1rm<- topTable(fit_m1mvl1rm, sort.by = "P", n = Inf)

#This is what the output column mean
#logFC: log2 fold change of I5.9/I5.6
#AveExpr: Average expression across all samples, in log2 CPM
#t: logFC divided by its standard error
#P.Value: Raw p-value (based on t) from test that logFC differs from 0
#adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
#B: log-odds that gene is DE (arguably less useful than the other columns)

#How many DE genes are there?
#lentil
length(which(toptab_l1lvl1rl$adj.P.Val < 0.05)) #189
length(which(toptab_l1lvl1m$adj.P.Val < 0.05)) #56
length(which(toptab_l1lvl1rm$adj.P.Val < 0.05)) #90
length(which(toptab_l1lvm1m$adj.P.Val < 0.05)) #1445

#mung
length(which(toptab_m1mvl1m$adj.P.Val < 0.05)) #285
length(which(toptab_m1mvl1rl$adj.P.Val < 0.05)) #1383
length(which(toptab_m1mvl1rm$adj.P.Val < 0.05)) #1579

#write out table to a file 
#Lentil
toptab_l1lvl1rl$Gene <- rownames(toptab_l1lvl1rl)
toptab_l1lvl1rl <- toptab_l1lvl1rl[,c("Gene", names(toptab_l1lvl1rl)[1:6])]
l1lvl1rl_deg<-toptab_l1lvl1rl[which(toptab_l1lvl1rl$adj.P.Val < 0.05),]
dim(l1lvl1rl_deg)
geneinfo1<-seqdata[which(l1lvl1rl_deg$Gene %in% seqdata$Geneid),][,c(2:4,30)]
l1lvl1rl_deg_inf<-cbind(l1lvl1rl_deg, geneinfo1)
head(l1lvl1rl_deg_inf)
dim(l1lvl1rl_deg_inf)
write.table(l1lvl1rl_deg_inf, file = "limma_results/L1LvL1RL_deg.txt", row.names = F, col.names = T, quote = F)

toptab_l1lvl1m$Gene <- rownames(toptab_l1lvl1m)
toptab_l1lvl1m <- toptab_l1lvl1m[,c("Gene", names(toptab_l1lvl1m)[1:6])]
l1lvl1m_deg<-toptab_l1lvl1m[which(toptab_l1lvl1m$adj.P.Val < 0.05),]
dim(l1lvl1m_deg)
geneinfo2<-seqdata[which(l1lvl1m_deg$Gene %in% seqdata$Geneid),][,c(2:4,30)]
l1lvl1m_deg_inf<-cbind(l1lvl1m_deg, geneinfo2)
head(l1lvl1m_deg_inf)
dim(l1lvl1m_deg_inf)
write.table(l1lvl1m_deg_inf, file = "limma_results/L1LvL1M_deg.txt", row.names = F, col.names = T, quote = F)

toptab_l1lvl1rm$Gene <- rownames(toptab_l1lvl1rm)
toptab_l1lvl1rm <- toptab_l1lvl1rm[,c("Gene", names(toptab_l1lvl1rm)[1:6])]
l1lvl1rm_deg<-toptab_l1lvl1rm[which(toptab_l1lvl1rm$adj.P.Val < 0.05),]
dim(l1lvl1rm_deg)
geneinfo3<-seqdata[which(l1lvl1rm_deg$Gene %in% seqdata$Geneid),][,c(2:4,30)]
l1lvl1rm_deg_inf<-cbind(l1lvl1rm_deg, geneinfo3)
dim(l1lvl1rm_deg_inf)
write.table(l1lvl1rm_deg_inf, file = "limma_results/L1LvL1RM_deg.txt", row.names = F, col.names = T, quote = F)

toptab_l1lvm1m$Gene <- rownames(toptab_l1lvm1m)
toptab_l1lvm1m <- toptab_l1lvm1m[,c("Gene", names(toptab_l1lvm1m)[1:6])]
l1lvm1m_deg<-toptab_l1lvm1m[which(toptab_l1lvm1m$adj.P.Val < 0.05),]
dim(l1lvm1m_deg)
geneinfo4<-seqdata[which(l1lvm1m_deg$Gene %in% seqdata$Geneid),][,c(2:4,30)]
l1lvm1m_deg_inf<-cbind(l1lvm1m_deg, geneinfo4)
dim(l1lvm1m_deg_inf)
write.table(l1lvm1m_deg_inf, file = "limma_results/L1LvM1M_deg.txt", row.names = F, col.names = T, quote = F)

#Mung
toptab_m1mvl1m$Gene <- rownames(toptab_m1mvl1m)
toptab_m1mvl1m <- toptab_m1mvl1m[,c("Gene", names(toptab_m1mvl1m)[1:6])]
m1mvl1m_deg<-toptab_m1mvl1m[which(toptab_m1mvl1m$adj.P.Val < 0.05),]
dim(m1mvl1m_deg)
geneinfo5<-seqdata[which(m1mvl1m_deg$Gene %in% seqdata$Geneid),][,c(2:4,30)]
m1mvl1m_deg_inf<-cbind(m1mvl1m_deg, geneinfo5)
dim(m1mvl1m_deg_inf)
write.table(m1mvl1m_deg_inf, file = "limma_results/M1MvL1M_deg.txt", row.names = F, col.names = T, quote = F)

toptab_m1mvl1rl$Gene <- rownames(toptab_m1mvl1rl)
toptab_m1mvl1rl <- toptab_m1mvl1rl[,c("Gene", names(toptab_m1mvl1rl)[1:6])]
m1mvl1rl_deg<-toptab_m1mvl1rl[which(toptab_m1mvl1rl$adj.P.Val < 0.05),]
dim(m1mvl1rl_deg)
geneinfo6<-seqdata[which(m1mvl1rl_deg$Gene %in% seqdata$Geneid),][,c(2:4,30)]
m1mvl1rl_deg_inf<-cbind(m1mvl1rl_deg, geneinfo6)
dim(m1mvl1rl_deg_inf)
write.table(m1mvl1rl_deg_inf, file = "limma_results/M1MvL1RL_deg.txt", row.names = F, col.names = T, quote = F)

toptab_m1mvl1rm$Gene <- rownames(toptab_m1mvl1rm)
toptab_m1mvl1rm <- toptab_m1mvl1rm[,c("Gene", names(toptab_m1mvl1rm)[1:6])]
m1mvl1rm_deg<-toptab_m1mvl1rm[which(toptab_m1mvl1rm$adj.P.Val < 0.05),]
dim(m1mvl1rm_deg)
geneinfo7<-seqdata[which(m1mvl1rm_deg$Gene %in% seqdata$Geneid),][,c(2:4,30)]
m1mvl1rm_deg_inf<-cbind(m1mvl1rm_deg, geneinfo7)
dim(m1mvl1rm_deg_inf)
write.table(m1mvl1rm_deg_inf, file = "limma_results/M1MvL1RM_deg.txt", row.names = F, col.names = T, quote = F)

##get summary of differentially expressed genes for each contrast
#lentil
summa_fit_l1lvl1rl<- decideTests(fit_l1lvl1rl)
summary(summa_fit_l1lvl1rl)

summa_fit_l1lvl1m<- decideTests(fit_l1lvl1m)
summary(summa_fit_l1lvl1m)

summa_fit_l1lvl1rm<- decideTests(fit_l1lvl1rm)
summary(summa_fit_l1lvl1rm)

summa_fit_l1lvm1m<- decideTests(fit_l1lvm1m)
summary(summa_fit_l1lvm1m)

#mung
summa_fit_m1mvl1m<- decideTests(fit_m1mvl1m)
summary(summa_fit_m1mvl1m)

summa_fit_m1mvl1rl<- decideTests(fit_m1mvl1rl)
summary(summa_fit_m1mvl1rl)

summa_fit_m1mvl1rm<- decideTests(fit_m1mvl1rm)
summary(summa_fit_m1mvl1rm)

#plots after testing for DE
#Let’s do a few plots to make sure everything looks good and that we haven’t made a mistake in the analysis. Genome-wide plots that are useful for checking are MAplots (or MDplots) and volcano plots. There are functions in limma for plotting these with fit.cont as input.

# We want to highlight the significant genes. We can get this from decideTests.
#Lentils
pdf("lentil_de_volcano_plot.pdf", width=10, height=10)
#L1LvL1RL
par(mfrow=c(4,2))
par(mar=c(5,5,2,4))
plotMD(fit_l1lvl1rl,coef=1,status=summa_fit_l1lvl1rl[,"groupL1.L - groupL1R.L"], values = c(-1, 1),cex=1, cex.lab=1, hl.col=brewer.pal(n = 3, name = "Dark2"))
# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 50 most DE genes
par(mar=c(5,5,2,4))
volcanoplot(fit_l1lvl1rl,coef=1,highlight=50,cex=1,cex.lab=1, names=row.names(coef(fit_l1lvl1rl)))

par(mar=c(5,5,2,4))
plotMD(fit_l1lvl1m,coef=1,status=summa_fit_l1lvl1m[,"groupL1.L - groupL1.M"], values = c(-1, 1),cex=1, cex.lab=1, hl.col=brewer.pal(n = 3, name = "Dark2"))
par(mar=c(5,5,2,4))
volcanoplot(fit_l1lvl1m,coef=1,highlight=50,cex=1,cex.lab=1, names=row.names(coef(fit_l1lvl1m)))

par(mar=c(5,5,2,4))
plotMD(fit_l1lvl1rm,coef=1,status=summa_fit_l1lvl1rm[,"groupL1.L - groupL1R.M"], values = c(-1, 1),cex=1, cex.lab=1, hl.col=brewer.pal(n = 3, name = "Dark2"))
par(mar=c(5,5,2,4))
volcanoplot(fit_l1lvl1rm,coef=1,highlight=50,cex=1,cex.lab=1, names=row.names(coef(fit_l1lvl1rm)))

par(mar=c(5,5,2,4))
plotMD(fit_l1lvm1m,coef=1,status=summa_fit_l1lvm1m[,"groupL1.L - groupM1.M"], values = c(-1, 1),cex=1, cex.lab=1, hl.col=brewer.pal(n = 3, name = "Dark2"))
par(mar=c(5,5,2,4))
volcanoplot(fit_l1lvm1m,coef=1,highlight=50,cex=1,cex.lab=1, names=row.names(coef(fit_l1lvm1m)))

dev.off()

#Mung
pdf("mung_de_volcano_plot.pdf", width=10, height=10)

par(mfrow=c(3,2))
par(mar=c(5,5,2,4))
plotMD(fit_m1mvl1m,coef=1,status=summa_fit_m1mvl1m[,"groupM1.M - groupL1.M"], values = c(-1, 1),cex=1, cex.lab=1, hl.col=brewer.pal(n = 3, name = "Dark2"))
par(mar=c(5,5,2,4))
volcanoplot(fit_m1mvl1m,coef=1,highlight=50,cex=1,cex.lab=1, names=row.names(coef(fit_m1mvl1m)))

plotMD(fit_m1mvl1rl,coef=1,status=summa_fit_m1mvl1rl[,"groupM1.M - groupL1R.L"], values = c(-1, 1),cex=1, cex.lab=1, hl.col=brewer.pal(n = 3, name = "Dark2"))
par(mar=c(5,5,2,4))
volcanoplot(fit_m1mvl1rl,coef=1,highlight=50,cex=1,cex.lab=1, names=row.names(coef(fit_m1mvl1rl)))

plotMD(fit_m1mvl1rm,coef=1,status=summa_fit_m1mvl1rm[,"groupM1.M - groupL1R.M"], values = c(-1, 1),cex=1, cex.lab=1, hl.col=brewer.pal(n = 3, name = "Dark2"))
par(mar=c(5,5,2,4))
volcanoplot(fit_m1mvl1rm,coef=1,highlight=50,cex=1,cex.lab=1, names=row.names(coef(fit_m1mvl1rm)))

dev.off()





