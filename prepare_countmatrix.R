"""
author: samridhi chaturvedi
samridhi.chaturvedi@gmail.com
data: March 2020
This script takes featurecount output files for each sample and combine the counts data of each sample into one count matrix.
"""

#read in all the files
samplefiles<-list.files(path=".", pattern = "Aligned.featureCounts$", all.files=T)

##order of samples 
[1] "cmac16Aligned.featureCounts"
[1] 26 78  0  0  0  0
[1] "cmac17Aligned.featureCounts"
[1]  2 32  0  0  0  0
[1] "cmac18Aligned.featureCounts"
[1] 34 48  0  0  0  0
[1] "cmac19Aligned.featureCounts"
[1] 22 30  0  0  0  0
[1] "cmac28Aligned.featureCounts"
[1] 32 46  0  0  0  0
[1] "cmac29Aligned.featureCounts"
[1] 30 38  0  0  0  0
[1] "cmac31Aligned.featureCounts"
[1] 12 54  2  0  0  0
[1] "cmac40Aligned.featureCounts"
[1] 10 72  0  0  0  0
[1] "cmac41Aligned.featureCounts"
[1] 20 40  0  0  0  0
[1] "cmac45Aligned.featureCounts"
[1] 22 82  0  0  0  0
[1] "cmac58Aligned.featureCounts"
[1] 44 70  0  0  0  0
[1] "cmac59Aligned.featureCounts"
[1] 24 56  0  0  0  0
[1] "cmac63Aligned.featureCounts"
[1] 38 20  1  0  0  0
[1] "cmac64Aligned.featureCounts"
[1] 16 62  0  0  0  0
[1] "cmac65Aligned.featureCounts"
[1] 48 44  0  0  0  0
[1] "cmac66Aligned.featureCounts"
[1] 50 32  0  0  0  0
[1] "cmac67Aligned.featureCounts"
[1] 36 40  0  0  0  0
[1] "cmac68Aligned.featureCounts"
[1] 30 30  0  0  0  0
[1] "cmac69Aligned.featureCounts"
[1] 54 88  0  0  0  0
[1] "cmac6Aligned.featureCounts"
[1] 18 34  1  0  0  0
[1] "cmac70Aligned.featureCounts"
[1] 22 50  0  0  0  0
[1] "cmac71Aligned.featureCounts"
[1] 30 20  0  0  0  0
[1] "cmac72Aligned.featureCounts"
[1] 28 48  0  0  0  0
[1] "cmac73Aligned.featureCounts"
[1] 36 52  1  0  0  0
[1] "cmac7Aligned.featureCounts"
[1] 42 64  0  0  0  0


#get dimensions of each file
for (file in samplefiles){
	sample<-read.table(file, skip=2)
	print (file)
	print (head(sample[,7]))
	}
	
#create the count matrix with all samples
genes <- read.table(samplefiles[1], skip=1, header=T)[,1:4]
df <- do.call(cbind,lapply(samplefiles,function(fn)read.table(fn,header=FALSE, skip=2, sep="\t")[,7]))
genes_df <- cbind(genes,df)
#name columns with sample ids
colnames(genes_df)[5:29]<-c("cmac16", "cmac17","cmac18","cmac19","cmac28","cmac29","cmac31","cmac40","cmac41","cmac45","cmac58","cmac59","cmac63","cmac64","cmac65","cmac66","cmac67","cmac68","cmac69","cmac6","cmac70","cmac71","cmac72","cmac73","cmac7")
#save file for further analyses
write.table(genes_df, "allsamples_featurecounts.txt", col.names=T, row.names=F, sep= " ", quote=F)

