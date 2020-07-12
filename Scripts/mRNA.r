###########################################
#  -Integrative bioinformatics            # 
#  - RNA-Seq                              #
#  - 2018- 12- 18                         #
#  - Copyright: Mohamed Hamed             #
###########################################
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# # Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

library(readr)
library(org.Hs.eg.db)
library("vsn")
library(DESeq2)



##### load the mRNA-Seq data #####
data.path="./Data/"
files <- list.files(path=data.path,recursive=T, pattern = ".gz")

# read the first file for the first time
file=files[1]
file.id=strsplit(file,"/")[[1]][1]

#open a connection to your gz file and read the file
gz.con=gzfile(file.path(data.path,files[1]))
temp <- read.table(gz.con, header=F)

#create a storing object mrna.exp to save the whole read counts of each file read in an iteration
mrna.exp=temp
rownames(mrna.exp)=mrna.exp[,1]
mrna.exp=mrna.exp[-1]
colnames(mrna.exp)=c(file.id)

for(i in 2: length(files))
{
  
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  
  # read the next file  
  gz.con=gzfile(file.path(data.path,files[i]))
  temp <- read.table(gz.con, header=F)
  
  ## remove the first column, bec we had it already
  temp=temp[-1]
  colnames(temp)=c(file.id)
  
  mrna.exp=cbind(mrna.exp,temp)
}
View(mrna.exp)

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

# check duplciation of of gene symbols?  
x=duplicated(mrna.exp2$symbol)  
sum(x)

### yes .. why ? transcripts?  solutions : aggregation
mrna.exp.data=mrna.exp2[-dim(mrna.exp2)[2]]
mrna.exp.data=apply(mrna.exp.data,2, as.numeric)

####remove  duplication by aggregation
mrna.exp.data.agg= aggregate(mrna.exp.data, list(mrna.exp2$symbol),FUN=mean)

rownames(mrna.exp.data.agg)=mrna.exp.data.agg$Group.1
mrna.exp.data.agg=mrna.exp.data.agg[-1]

file.ids=colnames(mrna.exp.data.agg)



###### load the mrna sample sheets


pheno <- read_delim("./Data/gdc_sample_sheet.2018-12-21.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
View(pheno)
table(pheno$`Sample Type`)

# rename column names : replace the spaces with dots 
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)


#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 

file.ids.pheno=pheno$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mrna.exp.data.agg)=pheno$Sample.ID[index.files]

phenob = pheno

# selecte only submitter ID and tumor stage then rename it to Case.ID then perform merging base on Case.ID
clinical_sub <- clinical$submitter_id
clinical_sub <- as.data.frame(clinical_sub)
clinical_sub <- cbind(clinical_sub, clinical$tumor_stage)
colnames(clinical_sub) = c("Case.ID", "Stage")
pheno$stage <- merge(pheno, clinical_sub, by= "Case.ID")
pheno$stage <- pheno$stage$Stage
pheno$Type.Stage <- paste(pheno$Sample.Type, pheno$stage)
pheno$Type.Stage = ifelse(pheno$Sample.Type == "Solid Tissue Normal", pheno$Type.Stage <- "Solid Tissue Normal", paste(pheno$Sample.Type, pheno$stage) )




# pheno=pheno[index.files,]
# names(mrna.exp.data.agg)= pheno$Sample.ID

# for simplifying the analysis (and for time considerations) we will consider only  20 sample from each type (normal and cancer)
# sample.no=20
all.normal.samples= pheno[ pheno$Type.Stage %in% c("Solid Tissue Normal"),]$Sample.ID
normal.samples=all.normal.samples

all.stageI.samples= pheno[ pheno$Type.Stage %in% c("Primary Tumor stage i"),]$Sample.ID
stageI.samples=all.stageI.samples

all.stageII.samples= pheno[ pheno$Type.Stage %in% c("Primary Tumor stage ii"),]$Sample.ID
stageII.samples=all.stageII.samples

all.stageIII.samples= pheno[ pheno$Type.Stage %in% c("Primary Tumor stage iii"),]$Sample.ID
stageIII.samples=all.stageIII.samples

all.stageIV.samples= pheno[ pheno$Type.Stage %in% c("Primary Tumor stage iv"),]$Sample.ID
stageIV.samples=all.stageIV.samples



# now we will retrieve  the exp data for these 20 samples only and also the pheno data
normal.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% normal.samples]
stageI.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% stageI.samples]
stageII.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% stageII.samples]
stageIII.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% stageIII.samples]
stageIV.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% stageIV.samples]

pheno.sub=pheno[pheno$Sample.ID %in% c(normal.samples,stageI.samples,stageII.samples,stageIII.samples,stageIV.samples), c("Sample.ID", "Type.Stage")]
#pheno.sub=pheno[pheno$Sample.ID %in% c(normal.samples,stageI.samples), c("Sample.ID", "Type.Stage")]




exp.sub=cbind(normal.exp,stageI.exp, stageII.exp, stageIII.exp, stageIV.exp)
exp.sub=apply (exp.sub, 2,as.integer)
rownames(exp.sub)=rownames(normal.exp)


save(mrna.exp.data.agg,pheno, exp.sub,pheno.sub ,file="RNA-seq.RDATA")

###### DO the differential EXP analysis using DeSeq2
cond1="Solid Tissue Normal" 
cond2="Primary Tumor stage i"
cond3="Primary Tumor stage ii"
cond4="Primary Tumor stage iii"
cond5="Primary Tumor stage iv"




dds = DESeqDataSetFromMatrix( countData = exp.sub , colData = pheno.sub , design = ~ Type.Stage)
dds.run = DESeq(dds)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res=results(dds.run)
res=results(dds.run, contrast = c("Type.Stage",cond1 ,cond2) )

# remove nulls
res=res[complete.cases(res), ]
summary(res)


res.df1_2=as.data.frame(res)

plotMA(res, ylim=c(-1,1)) 
summary (res)

res.degs1_2=res.df1_2[res.df1_2$padj< 0.05 & abs(res.df1_2$log2FoldChange)>log2(2),]


res2=results(dds.run, contrast = c("Type.Stage",cond2 ,cond3) )

# remove nulls
res2=res2[complete.cases(res2), ]
summary(res2)


res.df2_3=as.data.frame(res2)

plotMA(res2, ylim=c(-1,1)) 
summary (res2)

res.degs2_3=res.df2_3[res.df2_3$padj< 0.05 & abs(res.df2_3$log2FoldChange)>log2(2),]



res3=results(dds.run, contrast = c("Type.Stage",cond3 ,cond4) )

# remove nulls
res3=res3[complete.cases(res3), ]
summary(res3)


res.df3_4=as.data.frame(res3)

plotMA(res3, ylim=c(-1,1)) 
summary (res3)

res.degs3_4=res.df3_4[res.df3_4$padj< 0.05 & abs(res.df3_4$log2FoldChange)>log2(2),]


res4=results(dds.run, contrast = c("Type.Stage",cond4 ,cond5) )

# remove nulls
res4=res4[complete.cases(res4), ]
summary(res4)


res.df4_5=as.data.frame(res4)

plotMA(res4, ylim=c(-1,1)) 
summary (res4)

res.degs4_5=res.df4_5[res.df4_5$padj< 0.05 & abs(res.df4_5$log2FoldChange)>log2(2),]



#### get the normalized and loggedtransformed values of all exp data
#using the the variance stabilizing transformation. vsn package
ntd=normTransform(dds)
exp.sub.norm= assay(ntd)







# Make a basic volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


res.degs1_2$regulation = 1
res.degs1_2$regulation = ifelse(res.degs1_2$log2FoldChange > 0, 1,-1)
d1_2 = as.data.frame(res.degs1_2$regulation)
rownames(d1_2) = rownames(res.degs1_2)
write.table(d1_2, "mrna_DEGs1_2.txt", sep = "\t", col.names = FALSE)

res.degs2_3$regulation = 1
res.degs2_3$regulation = ifelse(res.degs2_3$log2FoldChange > 0, 1,-1)
d2_3 = as.data.frame(res.degs2_3$regulation)
rownames(d2_3) = rownames(res.degs2_3)
write.table(d2_3, "mrna_DEGs2_3.txt",  sep = "\t", col.names = FALSE)

res.degs3_4$regulation = 1
res.degs3_4$regulation = ifelse(res.degs3_4$log2FoldChange > 0, 1,-1)
d3_4 = as.data.frame(res.degs3_4$regulation)
rownames(d3_4) = rownames(res.degs3_4)
write.table(d3_4, "mrna_DEGs3_4.txt", sep = "\t", col.names = FALSE)

res.degs4_5$regulation = 1
res.degs4_5$regulation = ifelse(res.degs4_5$log2FoldChange > 0, 1,-1)
d4_5 = as.data.frame(res.degs4_5$regulation)
rownames(d4_5) = rownames(res.degs4_5)
write.table(d4_5, "mrna_DEGs4_5.txt",  sep = "\t", col.names = FALSE)


# 
# # label file ids of both tumor and normlal
# normal.samples.file.ids= pheno[ pheno$Sample.Type %in% c("Solid Tissue Normal"),]$File.ID
# tumor.file.ids= pheno[ pheno$Sample.Type %in% c("Primary Tumor"),]$File.ID
# 
# # get the tumor and normal dataframe
# mrna.exp.normal=mrna.exp.data.agg[, names(mrna.exp.data.agg) %in% normal.file.ids]
# mrna.exp.tumor=mrna.exp.data.agg[, names(mrna.exp.data.agg) %in% tumor.file.ids]
# 
# # rename column names of both tumor and normal
# tumor.names= paste (  rep("PCa",dim(mrna.exp.tumor)[2]) , 1:  dim(mrna.exp.tumor)[2], sep="-")
# normal.names= paste (  rep("Ctrl",dim(mrna.exp.normal)[2]) , 1:  dim(mrna.exp.normal)[2], sep="-")
# 
# sample.ids=pheno[pheno$File.ID %in% file.ids,]$Sample.ID
# names(mrna.exp.normal)=normal.names
# names(mrna.exp.tumor)=tumor.names



pheno$Type.Stage = ifelse(pheno$Sample.ID %in% c("TCGA-BL-A0C8-01A", "TCGA-BL-A13I-01A"), pheno$Type.Stage <- "noooo", ifelse(pheno$Sample.Type == "Solid Tissue Normal", pheno$Type.Stage <- "Solid Tissue Normal", paste(pheno$Sample.Type, pheno$stage) ) )


clinical <- read.delim("/Users/Roza/Documents/integrative_biology/bladder/DE/DEGs/clinical.tsv")
pheno = 
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)
