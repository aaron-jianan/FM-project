######Intersection genes
library(readr)
GSE221921_Sig <- read_csv("4.DEG_edgeR_Sig.csv")
GSE221921_Sig=na.omit(GSE221921_Sig)
DEGname=GSE221921_Sig$GeneID

library(readr)
GSE67311deg <- read_csv("GSE67311deg.csv")
GSE67311deg_SIG=GSE67311deg[GSE67311deg$P.Value<0.05,]
DEGname1=GSE67311deg_SIG$symbol

library(readr)
GSE229750fmhc_DEG_Sig <- read_csv("GSE229750fmhc_DEG_Sig.csv")
View(GSE229750fmhc_DEG_Sig)
DEGname2=GSE229750fmhc_DEG_Sig$SYMBOL
DEGname2=na.omit(DEGname2)

intername=intersect(intersect(DEGname,DEGname1),DEGname2)

##Venn
library (VennDiagram) 
venn.diagram(x=list(DEGname,DEGname1,DEGname2),
             scaled = F, 
             alpha= 0.5, 
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), 
             label.col ='black' ,
             cex = 2, 
             fontface = "bold",  
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), 
             category.names = c("GSE221921", "GSE67311","GSE229750") , 
             cat.dist = 0, 
             cat.pos = c(-120, -240, -180), 
             cat.cex = 1, 
             cat.fontface = "bold", 
             cat.col='black' ,   
             cat.default.pos = "outer",
             output=TRUE,
             filename='venn.png',
             imagetype="png", 
             resolution = 400,  
             compression = "lzw")

#########################GSE221921###########################################
load("20240701GSE221921.Rdata")
rm(list = setdiff(ls(), c("filter_count")))
filter_count=as.data.frame(t(filter_count))
genes_of_interest <- c("ARHGEF37", "DYRK3", "RGS17")

##COR>0.6 P<0.05
cor1=cor(filter_count,method = c("spearman"))
cor2=cor1[,colnames(cor1) %in% genes_of_interest]

RGS17_GSE221921=as.data.frame(cor2[,1])
RGS17_GSE221921$name=rownames(RGS17_GSE221921)
DYRK3_GSE221921=as.data.frame(cor2[,2])
DYRK3_GSE221921$name=rownames(DYRK3_GSE221921)
ARHGEF37_GSE221921=as.data.frame(cor2[,3])
ARHGEF37_GSE221921$name=rownames(ARHGEF37_GSE221921)


RGS17_GSE221921_0.6=as.data.frame(RGS17_GSE221921[abs(RGS17_GSE221921$`cor2[, 1]`)>0.6,])
RGS17_name=rownames(RGS17_GSE221921_0.6)

DYRK3_GSE221921_0.6=as.data.frame(DYRK3_GSE221921[abs(DYRK3_GSE221921$`cor2[, 2]`)>0.6,])
DYRK3_name=rownames(DYRK3_GSE221921_0.6)

ARHGEF37_GSE221921_0.6=as.data.frame(ARHGEF37_GSE221921[abs(ARHGEF37_GSE221921$`cor2[, 3]`)>0.6,])
ARHGEF37_name=rownames(ARHGEF37_GSE221921_0.6)
RGS17_COR=filter_count[,colnames(filter_count) %in% RGS17_name]
DYRK3_COR=filter_count[,colnames(filter_count) %in% DYRK3_name]
ARHGEF37_COR=filter_count[,colnames(filter_count) %in% ARHGEF37_name]

library(Hmisc)
RGS17=rcorr(as.matrix(RGS17_COR), type = "spearman")
View(RGS17$r)
RGS17$r=round(RGS17$r,2)
View(RGS17$P)

DYRK3=rcorr(as.matrix(DYRK3_COR), type = "spearman")
View(DYRK3$r)
DYRK3$r=round(DYRK3$r,2)
View(DYRK3$P)

ARHGEF37=rcorr(as.matrix(ARHGEF37_COR), type = "spearman")
View(ARHGEF37$r)
ARHGEF37$r=round(ARHGEF37$r,2)
View(ARHGEF37$P)
ARHGEF37 <-as.data.frame(ARHGEF37$P)
colnames(ARHGEF37)
ARHGEF37=ARHGEF37[ARHGEF37$ARHGEF37<0.05,]
rownames(ARHGEF37)
ARHGEF37_name=rownames(ARHGEF37)

DYRK3 <-as.data.frame(DYRK3$P)
colnames(DYRK3)
DYRK3=DYRK3[DYRK3$DYRK3<0.05,]
rownames(DYRK3)
DYRK3_name=rownames(DYRK3)

RGS17 <-as.data.frame(RGS17$P)
colnames(RGS17)
RGS17=RGS17[RGS17$RGS17<0.05,]
rownames(RGS17)
RGS17_name=rownames(RGS17)

RGS17_name ##343
DYRK3_name ##351
ARHGEF37_name #26

save(RGS17_name,DYRK3_name,ARHGEF37_name,file="GSE221921_genes.Rdata")

#########################GSE67311###########################################
rm(list = ls())
options(stringsAsFactors = F)
load(file="step1output.Rdata")
library(stringr)
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL11532')
exp=as.data.frame(exp)
exp$symbol=rownames(exp)
exp1=merge(exp,ids,by.x="symbol",by.y="probe_id",all.x = T)
exp1=exp1[!is.na(exp1$symbol.y),]
exp1=exp1[!duplicated(exp1$symbol.y),]
rownames(exp1)=exp1$symbol.y
colnames(exp1)
exp1=exp1[,-c(1,144)]
exp1=as.data.frame(t(exp1))

rm(list = setdiff(ls(), c("exp1")))
genes_of_interest <- c("ARHGEF37", "DYRK3", "RGS17")
filter_count=exp1

##COR>0.6
cor1=cor(filter_count,method = c("spearman"))
cor2=cor1[,colnames(cor1) %in% genes_of_interest]

DYRK3_GSE67311=as.data.frame(cor2[,1])
DYRK3_GSE67311$name=rownames(DYRK3_GSE67311)
ARHGEF37_GSE67311=as.data.frame(cor2[,2])
ARHGEF37_GSE67311$name=rownames(ARHGEF37_GSE67311)

DYRK3_GSE67311_0.6=as.data.frame(DYRK3_GSE67311[abs(DYRK3_GSE67311$`cor2[, 1]`)>0.6,])
DYRK3_name=rownames(DYRK3_GSE67311_0.6)

ARHGEF37_GSE67311_0.6=as.data.frame(ARHGEF37_GSE67311[abs(ARHGEF37_GSE67311$`cor2[, 2]`)>0.6,])
ARHGEF37_name=rownames(ARHGEF37_GSE67311_0.6)

DYRK3_COR=filter_count[,colnames(filter_count) %in% DYRK3_name]
ARHGEF37_COR=filter_count[,colnames(filter_count) %in% ARHGEF37_name]

library(Hmisc)
DYRK3=rcorr(as.matrix(DYRK3_COR), type = "spearman")
View(DYRK3$r)
DYRK3$r=round(DYRK3$r,2)
View(DYRK3$P)

ARHGEF37=rcorr(as.matrix(ARHGEF37_COR), type = "spearman")
View(ARHGEF37$r)
ARHGEF37$r=round(ARHGEF37$r,2)
View(ARHGEF37$P)
ARHGEF37 <-as.data.frame(ARHGEF37$P)
colnames(ARHGEF37)
ARHGEF37=ARHGEF37[ARHGEF37$ARHGEF37<0.05,]
rownames(ARHGEF37)
ARHGEF37_name1=rownames(ARHGEF37)

DYRK3 <-as.data.frame(DYRK3$P)
colnames(DYRK3)
DYRK3=DYRK3[DYRK3$DYRK3<0.05,]
rownames(DYRK3)
DYRK3_name1=rownames(DYRK3)
DYRK3_name1=DYRK3_name ##218
ARHGEF37_name1=ARHGEF37_name #137

save(DYRK3_name1,ARHGEF37_name1,file="GSE67311_genes.Rdata")

####################GSE229750###################
load("20240701GSE229750.Rdata")
rm(list = setdiff(ls(), c("filter_count")))
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(filter_count), 
                  fromType = "ENTREZID", 
                  toType = "SYMBOL", 
                  OrgDb = org.Hs.eg.db)
head(id2symbol)
filter_count$ENTREZID=rownames(filter_count)
filter_count <- merge(id2symbol,filter_count,
                          by.x="ENTREZID",by.y="ENTREZID",all.y=T)

filter_count=filter_count[!is.na(filter_count$SYMBOL),]
filter_count=filter_count[!duplicated(filter_count$SYMBOL),]
rownames(filter_count)=filter_count$SYMBOL
colnames(filter_count)
filter_count=filter_count[,-c(1,2)]
filter_count=filter_count[,-c(2,4)]
filter_count=as.data.frame(t(filter_count))
genes_of_interest <- c("ARHGEF37", "DYRK3", "RGS17")

cor1=cor(filter_count,method = c("spearman"))
cor2=cor1[,colnames(cor1) %in% genes_of_interest]

RGS17_GSE229750=as.data.frame(cor2[,1])
RGS17_GSE229750$name=rownames(RGS17_GSE229750)
DYRK3_GSE229750=as.data.frame(cor2[,3])
DYRK3_GSE229750$name=rownames(DYRK3_GSE229750)
ARHGEF37_GSE229750=as.data.frame(cor2[,2])
ARHGEF37_GSE229750$name=rownames(ARHGEF37_GSE229750)


RGS17_GSE229750_0.6=as.data.frame(RGS17_GSE229750[abs(RGS17_GSE229750$`cor2[, 1]`)>0.6,])
RGS17_name=rownames(RGS17_GSE229750_0.6)

DYRK3_GSE229750_0.6=as.data.frame(DYRK3_GSE229750[abs(DYRK3_GSE229750$`cor2[, 3]`)>0.6,])
DYRK3_name=rownames(DYRK3_GSE229750_0.6)

ARHGEF37_GSE229750_0.6=as.data.frame(ARHGEF37_GSE229750[abs(ARHGEF37_GSE229750$`cor2[, 2]`)>0.6,])
ARHGEF37_name=rownames(ARHGEF37_GSE229750_0.6)

RGS17_COR=filter_count[,colnames(filter_count) %in% RGS17_name]
DYRK3_COR=filter_count[,colnames(filter_count) %in% DYRK3_name]
ARHGEF37_COR=filter_count[,colnames(filter_count) %in% ARHGEF37_name]

library(Hmisc)
RGS17=rcorr(as.matrix(RGS17_COR), type = "spearman")
View(RGS17$r)
RGS17$r=round(RGS17$r,2)
View(RGS17$P)

DYRK3=rcorr(as.matrix(DYRK3_COR), type = "spearman")
View(DYRK3$r)
DYRK3$r=round(DYRK3$r,2)
View(DYRK3$P)

ARHGEF37=rcorr(as.matrix(ARHGEF37_COR), type = "spearman")
View(ARHGEF37$r)
ARHGEF37$r=round(ARHGEF37$r,2)
View(ARHGEF37$P)

ARHGEF37 <-as.data.frame(ARHGEF37$P)
colnames(ARHGEF37)
ARHGEF37=ARHGEF37[ARHGEF37$ARHGEF37<0.05,]
ROWNAMES(ARHGEF37)
ARHGEF37_names2=ROWNAMES(ARHGEF37)

DYRK3 <-as.data.frame(DYRK3$P)
colnames(DYRK3)
DYRK3=DYRK3[DYRK3$DYRK3<0.05,]
ROWNAMES(DYRK3)
DYRK3_names2=ROWNAMES(DYRK3)

RGS17 <-as.data.frame(RGS17$P)
colnames(RGS17)
RGS17=RGS17[RGS17$RGS17<0.05,]
ROWNAMES(RGS17)
RGS17_names2=ROWNAMES(RGS17)

save(RGS17_names2,DYRK3_names2,ARHGEF37_names2,file="GSE229750_genes.Rdata")


#intersect
rm(list = ls())
load("GSE229750_genes.Rdata")
load("GSE67311_genes.Rdata")
load("GSE221921_genes.Rdata")

##Merge and remove duplicates
ARHGEF37_IMGENE=c(ARHGEF37_names2,ARHGEF37_name,ARHGEF37_name1)
ARHGEF37_IMGENE=na.omit(ARHGEF37_IMGENE)
ARHGEF37_IMGENE=unique(ARHGEF37_IMGENE)

DYRK3_IMGENE=c(DYRK3_names2,DYRK3_name,DYRK3_name1)
DYRK3_IMGENE=na.omit(DYRK3_IMGENE)
DYRK3_IMGENE=unique(DYRK3_IMGENE)

RGS17_IMGENE=c(RGS17_names2,RGS17_name)
RGS17_IMGENE=na.omit(RGS17_IMGENE)
RGS17_IMGENE=unique(RGS17_IMGENE)

save(ARHGEF37_IMGENE,DYRK3_IMGENE,RGS17_IMGENE,
     file = "IMgene.Rdata")

write.csv(ARHGEF37_IMGENE, file = "ARHGEF37_IMGENE.csv")
write.csv(DYRK3_IMGENE, file = "DYRK3_IMGENE.csv")
write.csv(RGS17_IMGENE, file = "RGS17_IMGENE.csv")