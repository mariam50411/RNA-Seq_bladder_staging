#load needed libraries
library(readr)
library(org.Hs.eg.db)
library("vsn")
library(plyr)
library(NMF)
library(gplots)
library(GOplot)
library(ggplot2)
library(DESeq2)

#load mRNA counts data
data.path= "E:/Usb/Downloads/Downloads/Bladder_stagess/RNA"
files <- list.files(path=data.path,recursive=T, pattern = "gz")
#make quantification matrix from all samples
file=files[1]
file.id=strsplit(file,"/")[[1]][1]
gz.con=gzfile(file.path(data.path,files[1]))
temp=read.table(gz.con, header=F)
mrna.exp=temp
rownames(mrna.exp)=mrna.exp[,1]
mrna.exp=mrna.exp[-1]
colnames(mrna.exp)=c(file.id)
for(i in 2: length(files))
{
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  gz.con=gzfile(file.path(data.path,files[i]))
  temp <- read.table(gz.con, header=F)
  temp=temp[-1]
  colnames(temp)=c(file.id)
  
  mrna.exp=cbind(mrna.exp,temp)
}
#View(mrna.exp)

### do the mapping of ensembel.id to gene symbol###############

# prepare the ensembel id to be like the one in the database
ensemble.id=sapply(rownames(mrna.exp), function(x) strsplit(as.character(x),"\\.")[[1]][1])
mrna.exp=cbind(ensemble.id,mrna.exp)

mapper<- mapIds(org.Hs.eg.db, keys=ensemble.id, column="SYMBOL",keytype="ENSEMBL", multiVals="first")
mapper.df=as.data.frame(mapper)

mapper.df=cbind(rownames(mapper.df), mapper.df)
names(mapper.df)=c("ensemble.id","symbol")

mrna.exp2=merge(mrna.exp,mapper.df,by="ensemble.id",all.x=T)
# drop the first column (ensemble.id)
rownames(mrna.exp2) = mrna.exp2[,1]
mrna.exp2=mrna.exp2[-1]

mrna.exp2=mrna.exp2[ ! is.na(mrna.exp2$symbol),]
sum(duplicated(mrna.exp2$symbol))

### yes .. why ? transcripts?  solutions : aggregation
mrna.exp.data=mrna.exp2[-dim(mrna.exp2)[2]]
mrna.exp.data=apply(mrna.exp.data,2, as.numeric)

####remove  duplication by aggregation
mrna.exp.data.agg= aggregate(mrna.exp.data, list(mrna.exp2$symbol),FUN=mean)

rownames(mrna.exp.data.agg)=mrna.exp.data.agg$Group.1
mrna.exp.data.agg=mrna.exp.data.agg[-1]

file.ids=colnames(mrna.exp.data.agg)


#### load the mrna sample sheets ####
pheno <- read_delim("gdc_sample_sheet.2020-07-13.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

# rename column names : replace the spaces with dots 
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)

#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 
file.ids.pheno=pheno$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mrna.exp.data.agg)=pheno$Sample.ID[index.files]


#### load the mrna clinical sheets which contains info about the stages##
clinical <- read.delim("E:/Usb/Downloads/Downloads/Bladder_stagess/RNA/clinical.tsv")
#remove duplicated cases rom the clinical sheet
clinical = clinical[!duplicated(clinical$case_submitter_id),]
# selecte only submitter ID and tumor stage then rename it 
#to Case.ID then perform merging base on Case.ID
clinical<- clinical[,c("case_submitter_id","tumor_stage")]
colnames(clinical) = c("Case.ID", "Stages")

#subset the sample sheet to normal samples(no info in the clinical sheet)
#and tumor samples (present in the clinical sheet with each stage) so
#we can merge betweent the clinical file and the tumor samples
pheno.tum=as.data.frame(pheno[pheno$Sample.Type == "Primary Tumor",])
pheno.norm=pheno[pheno$Sample.Type == "Solid Tissue Normal",]
#remove duplicated cases from tumor samples sheet
pheno.tum=as.data.frame(pheno.tum[!duplicated(pheno.tum$Case.ID),])
#merge between the clinical file and the tumor samples sheet
pheno.tum <- join(pheno.tum, clinical, type = "inner")
#put the the Sample Type and its stage in one column
pheno.tum$Type.Stage <- paste(pheno.tum$Sample.Type, pheno.tum$Stages)
#remove the coulmns of Sample.Type and stages
pheno.tum=pheno.tum[,c(-8,-9)]
#rename the Type.Stage coulmn into Sample.Type so we can bind it to
#the normal samples sheet
colnames(pheno.tum)[8]="Sample.Type"
#bind the new tumor sample sheet with contains info about stages to 
#the normal sample sheet
pheno=rbind(pheno.norm,pheno.tum)

#export the expression matrix and the pheno files 
write.csv(as.data.frame(mrna.exp.data.agg),file="exp.csv", quote=F,row.names=T)
write.csv(as.data.frame(pheno),file="pheno.csv", quote=F,row.names=T)
#clear the environment
##################################################################################
#load the data
exp = read.csv("E:/Usb/Downloads/Downloads/Bladder_stagess/RNA/exp.csv", row.names=1)
pheno = read.csv("E:/Usb/Downloads/Downloads/Bladder_stagess/RNA/pheno.csv")

#make the sample IDs row names of the pheno sheet
row.names(pheno)=pheno$Sample.ID
#subset the pheno sheet for each group 
norm = pheno[pheno$Sample.Type == "Solid Tissue Normal",]
sp1 = pheno[pheno$Sample.Type == "Primary Tumor stage i",]
sp2 = pheno[pheno$Sample.Type == "Primary Tumor stage ii",]
sp3 = pheno[pheno$Sample.Type == "Primary Tumor stage iii",]
sp4 = pheno[pheno$Sample.Type == "Primary Tumor stage iv",]

#subset the expression matrix for each group 
exp.norm = exp[ ,which((names(exp) %in% norm$Sample.ID) == TRUE)]
exp.sp1 = exp[ ,which((names(exp) %in% sp1$Sample.ID) ==TRUE)]
exp.sp2 = exp[ ,which((names(exp) %in% sp2$Sample.ID) ==TRUE)]
exp.sp3 = exp[ ,which((names(exp) %in% sp3$Sample.ID) ==TRUE)]
exp.sp4 = exp[ ,which((names(exp) %in% sp4$Sample.ID) ==TRUE)]
#make sure that both expression and pheno files contains only data
#with stages
pheno = rbind(norm,sp1,sp2,sp3,sp4)
exp= cbind(exp.norm,exp.sp1,exp.sp2,exp.sp3,exp.sp4)

#make sure the the row names of the pheno in the same order for the
#column names of the expression matrix
rownames(pheno)=pheno$Sample.ID
pheno=pheno[colnames(exp),]

#convert the data values to integers
genes=row.names(exp)
exp=apply(exp,2,as.integer)
row.names(exp)=genes

###### DO the differential EXP analysis using DeSeq2

#specify how many conditions do you want to compare according to 
#the phenotypic table
cond1="Solid Tissue Normal" 
cond2="Primary Tumor stage i"
cond3="Primary Tumor stage ii"
cond4="Primary Tumor stage iii"
cond5="Primary Tumor stage iv"

#creat a deseq dataset object
dds = DESeqDataSetFromMatrix( countData = exp , colData = pheno , design = ~ Sample.Type)
#run the deseq2 worflow
dds.run = DESeq(dds)
#specifying teh contrast (to make a res object based on two specific conditions)
res1=results(dds.run, contrast = c("Sample.Type",cond2 ,cond1) )
res2=results(dds.run, contrast = c("Sample.Type",cond3 ,cond2) )
res3=results(dds.run, contrast = c("Sample.Type",cond4 ,cond3) )
res4=results(dds.run, contrast = c("Sample.Type",cond5 ,cond4) )
res5=results(dds.run, contrast = c("Sample.Type",cond3 ,cond1) )

# remove nulls
res1=as.data.frame(res1[complete.cases(res1), ])
res2=as.data.frame(res2[complete.cases(res2), ])
res3=as.data.frame(res3[complete.cases(res3), ])
res4=as.data.frame(res4[complete.cases(res4), ])
res5=as.data.frame(res5[complete.cases(res5), ])

##chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 4
dg0.1=res1[res1$padj< 0.05 & abs(res1$log2FoldChange)>log2(4),]
dg0.2=res5[res5$padj< 0.05 & abs(res5$log2FoldChange)>log2(4),]
dg1.2=res2[res2$padj< 0.05 & abs(res2$log2FoldChange)>log2(4),]
dg2.3=res3[res3$padj< 0.05 & abs(res3$log2FoldChange)>log2(4),]
dg3.4=res4[res4$padj< 0.05 & abs(res4$log2FoldChange)>log2(4),]

#add regulation coulmn(1 for upregulated genes and -1 for down regulated genes)
dg0.1$regulation = 1
dg0.1$regulation = ifelse(dg0.1$log2FoldChange > 0, 1,-1)

dg0.2$regulation = 1
dg0.2$regulation = ifelse(dg0.2$log2FoldChange > 0, 1,-1)

dg1.2$regulation = 1
dg1.2$regulation = ifelse(dg1.2$log2FoldChange > 0, 1,-1)

dg2.3$regulation = 1
dg2.3$regulation = ifelse(dg2.3$log2FoldChange > 0, 1,-1)

dg3.4$regulation = 1
dg3.4$regulation = ifelse(dg3.4$log2FoldChange > 0, 1,-1)

#export the Degs into your current folder for further analysthis
write.csv(as.data.frame(dg0.1),file="dg0.1.csv", quote=F,row.names=T)
write.csv(as.data.frame(dg0.2),file="dg0.2.csv", quote=F,row.names=T)
write.csv(as.data.frame(dg1.2),file="dg1.2.csv", quote=F,row.names=T)
write.csv(as.data.frame(dg2.3),file="dg2.3.csv", quote=F,row.names=T)
write.csv(as.data.frame(dg3.4),file="dg3.4.csv", quote=F,row.names=T)

#drow DEGs volcano plot
par(mfrow=c(2,2))
with(res1, plot(log2FoldChange, -log10(padj), pch=20, main="Control vs Stage1 DEGs"))#, xlim=c(-10,11), ylim=c(0,10)))
with(subset(res1, padj<.05 & (log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res1, padj<.05 & (log2FoldChange)< -2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
#with(subset(result, pvalue< 0.05 & abs(log2FoldChange)>1.2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
legend(x=-20,y=30,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)


with(res2, plot(log2FoldChange, -log10(padj), pch=20, main="Stage1 vs Stage2 DEGs"))#, xlim=c(-8,11), ylim=c(0,8)))
with(subset(res2, padj<.05 & (log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res2, padj<.05 & (log2FoldChange)< -2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
#with(subset(result, pvalue< 0.05 & abs(log2FoldChange)>1.2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
legend(x=16.5,y=23,c("upregulated","downgulated"),cex=.8,col=c("blue","red"), bty="n",pch=19)


# with(res5, plot(log2FoldChange, -log10(padj), pch=20, main="Control vs Stage2 DEGs"))#, xlim=c(-10,11), ylim=c(0,10)))
# with(subset(res5, padj<.05 & (log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
# with(subset(res5, padj<.05 & (log2FoldChange)< -2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
# #with(subset(result, pvalue< 0.05 & abs(log2FoldChange)>1.2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# legend(x=25,y=165,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)

with(res3, plot(log2FoldChange, -log10(padj), pch=20, main="Stage2 vs Stage3 DEGs"))#, xlim=c(-4.7,5.7), ylim=c(0,9)))
with(subset(res3, padj<.05 & (log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res3, padj<.05 & (log2FoldChange)< -2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
#with(subset(result, pvalue< 0.05 & abs(log2FoldChange)>1.2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
legend(x=-3.5,y=23,c("upregulated","downgulated"),cex=.8, bty="n", col=c("blue","red"),pch=19)

with(res4, plot(log2FoldChange, -log10(padj), pch=20, main="Stage3 vs Stage4 DEGs"))#, xlim=c(-5,4.5), ylim=c(0,12)))
with(subset(res4, padj<.05 & (log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res4, padj<.05 & (log2FoldChange)< -2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
#with(subset(result, pvalue< 0.05 & abs(log2FoldChange)>1.2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
legend(x=1,y=17,c("upregulated","downgulated"),cex=.8, bty="n",col=c("blue","red"),pch=19)

####drow heatmap####
#normalize the data
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))

#subset the normalized counts matrix for each group 
exp.norm = normalized_counts[ ,which((names(normalized_counts) %in% norm$Sample.ID) == TRUE)]
exp.sp1 = normalized_counts[ ,which((names(normalized_counts) %in% sp1$Sample.ID) ==TRUE)]
exp.sp2 = normalized_counts[ ,which((names(normalized_counts) %in% sp2$Sample.ID) ==TRUE)]
exp.sp3 = normalized_counts[ ,which((names(normalized_counts) %in% sp3$Sample.ID) ==TRUE)]
exp.sp4 = normalized_counts[ ,which((names(normalized_counts) %in% sp4$Sample.ID) ==TRUE)]

#bind the pairwise comparesions groups for both count and pheno data
exp1=cbind(exp.norm,exp.sp1)
exp2=cbind(exp.sp1,exp.sp2)
exp3=cbind(exp.sp2,exp.sp3)
exp4=cbind(exp.sp3,exp.sp4)
exp5=cbind(exp.norm,exp.sp2)
#pheno groups
sp.1=rbind(norm,sp1)
sp.2=rbind(sp1,sp2)
sp.3=rbind(sp2,sp3)
sp.4=rbind(sp3,sp4)
sp.5=rbind(norm,sp2)

#extract counts values of DEGs only for each stage
exp.degs1=as.matrix(exp1[rownames(exp1) %in% rownames(dg0.1), ])
exp.degs5=as.matrix(exp5[rownames(exp5) %in% rownames(dg0.2), ])
exp.degs2=as.matrix(exp2[rownames(exp2) %in% rownames(dg1.2), ])
exp.degs3=as.matrix(exp3[rownames(exp3) %in% rownames(dg2.3), ])
exp.degs4=as.matrix(exp4[rownames(exp4) %in% rownames(dg3.4), ])

#drow the heatmap
par(mfrow=c(2,2))
aheatmap(log2(exp.degs1+1), annCol =sp.1$Sample.Type, col = rev(brewer.pal(9,"RdBu")), main="mRNA Control vs Stage1")
aheatmap(log2(exp.degs2+1),annCol =sp.2$Sample.Type,col = rev(brewer.pal(9,"RdBu")),main="mRNA Stage1 vs Stage2")
#aheatmap(log2(exp.degs5+1),annCol =sp.5$Sample.Type,col = rev(brewer.pal(9,"RdBu")), main="mRNA Control vs Stage2")
aheatmap(log2(exp.degs3+1),annCol =sp.3$Sample.Type,col = rev(brewer.pal(9,"RdBu")),main="mRNA Stage2 vs Stage3")
aheatmap(log2(exp.degs4+1),annCol =sp.4$Sample.Type,col = rev(brewer.pal(9,"RdBu")),main="mRNA Stage3 vs Stage4")


