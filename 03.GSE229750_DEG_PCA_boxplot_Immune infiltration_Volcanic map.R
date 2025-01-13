##Raw data processing
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(data.table)

counts<- fread("GSE229750_raw_counts_GRCh38.p13_NCBI.tsv")
library(GEOquery)
gse_number = "GSE229750"
eSet <- getGEO(gse_number, destdir = '.', getGPL = F)
pd <- pData(eSet)
p = identical(rownames(pd),colnames(counts)[2:13]);p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
gpl_number <- eSet@annotation;gpl_number

colnames( pd)
group <- pd[, c("geo_accession","title")]
group

rawcount=counts
rawcount=as.data.frame(rawcount)
rownames(rawcount)=counts$GeneID
rawcount=rawcount[,-1]

keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
table(keep)
filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)
library(edgeR)
express_cpm <- cpm(filter_count)
express_cpm[1:6,1:6]

library(FactoMineR)
library(factoextra)
library(corrplot)
library(pheatmap)
library(tidyverse)

# PCA
library(tidyverse)
group=group[-c(2,4),]
group$group=c("FM",
              "FM",
              "FM","FM",
              "FM","Control","Control","Control","Control","Control")

group$group=as.factor(group$group)
dat <- log10(express_cpm+1)
dat=dat[,-c(2,4)]
dat[1:4,1:4]
dim(dat)
sampleTree <- hclust(dist(t(dat)), method = "average")
plot(sampleTree)

temp <- as.data.frame(cutree(sampleTree,k = 2)) %>% 
  rownames_to_column(var="sample")
temp1 <- merge(temp,group,by.x = "sample",by.y="geo_accession")
table(temp1$`cutree(sampleTree, k = 2)`,temp1$group)

dat <- log10(express_cpm+1)
dat[1:4,1:4]
dat <- as.data.frame(t(dat))
dat_pca <- PCA(dat, graph = FALSE)
dat[1:4,1:4]

group_list <- group[match(group$geo_accession,rownames(dat)),3]
group_list
mythe <- theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

p <- fviz_pca_ind(dat_pca,
                  geom.ind = "text", 
                  col.ind = group_list, 
                  palette = c("#00AFBB", "#E7B800","#114454"), 
                  addEllipses = T,  
                  legend.title = "Groups") + mythe
p

# DEG
library(edgeR)
library(ggplot2)

filter_count[1:4,1:4]
library(tidyverse)
group$group=c("FMW0","FMW12",
              "FMW0","FMW12",
              "FMW0","FMW0",
              "FMW0","Control","Control","Control","Control","Control")

group_list=group[c(1,3,5:12),]
group_list$group=as.factor(group_list$group)
filter_count=filter_count[,-c(2,4)]

group_list <- group_list[match(colnames(filter_count),group_list$geo_accession),3]
table(group_list)

design <- model.matrix(~0+group_list)
rownames(design) <- colnames(filter_count)
colnames(design) <- levels(factor(group_list))
design

DEG <- DGEList(counts=filter_count, 
               group=factor(group_list))
DEG <- calcNormFactors(DEG)
DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)
fit <- glmFit(DEG, design)
lrt <- glmLRT(fit, contrast=c(-1,1)) 

DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
head(DEG_edgeR)

fc_cutoff <- 1.5
pvalue <- 0.05

DEG_edgeR$regulated <- "normal"

loc_up <- intersect(which( DEG_edgeR$logFC > log2(fc_cutoff) ),
                    which( DEG_edgeR$PValue < pvalue) )

loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))

DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_edgeR), 
                  fromType = "ENTREZID", 
                  toType = "SYMBOL", 
                  OrgDb = org.Hs.eg.db)
head(id2symbol)
DEG_edgeR=DEG_edgeR[,2:7]
colnames(DEG_edgeR)
DEG_edgeR$ENTREZID=rownames(DEG_edgeR)
DEG_edgeR_symbol <- merge(id2symbol,DEG_edgeR,
                          by.x="ENTREZID",by.y="ENTREZID",all.y=T)
head(DEG_edgeR_symbol)
library(tidyverse)
DEG_edgeR_Sig <- filter(DEG_edgeR_symbol,regulated!="normal")

write.csv(DEG_edgeR,"GSE229750fmhc_DEG_all.csv", row.names = F)
write.csv(DEG_edgeR_Sig,"GSE229750fmhc_DEG_Sig.csv", row.names = F)

##Volcanic map

library(ggplot2)
library(tidyverse)
data <- DEG_edgeR_symbol
colnames(data)
p <- ggplot(data=data, aes(x=logFC, y=-log10(PValue),color=regulated)) + 
  geom_point(alpha=0.5, size=1.2) + 
  theme_set(theme_set(theme_bw(base_size=20))) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c(down='blue',normal='grey',up='red')) +
  geom_vline(xintercept=c(-(log2(1.5)),log2(1.5)),lty=2,col="black",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.6)
p

##DYRK3/RGS17/ARHGEF37
which(data$SYMBOL=="DYRK3")#16869
which(data$SYMBOL=="RGS17")#8110
which(data$SYMBOL=="ARHGEF37")#10069

label <- data[c(16869,8110,10069),]
p1 <- p + geom_point(size = 1, shape = 2, data = label,color="red") +
  ggrepel::geom_text_repel( aes(label = SYMBOL), data = label, color="black" )

p1

png(file = "5.GSE2297550Volcano_Plot.png",width = 900, height = 800, res=150)
plot(p1)
dev.off()

##Expression box diagram
###8444  DYRK3
#	26575 RGS17
#389337 ARHGEF37

which(rownames(express_cpm)==8444) ##1713
which(rownames(express_cpm)==26575) ##6516
which(rownames(express_cpm)==389337) ##5417

exp1=express_cpm[c(1713,6516,5417),]
rownames(exp1)=c("DYRK3","RGS17","ARHGEF37")
exp1=as.data.frame(t(exp1))
exp1$group=group$group
exp1=exp1[-c(2,4),]
library(tidyverse)
exp1$group=str_replace_all(exp1$group,"FMW0","FM")

library(reshape)
library(ggplot2)
exp1=melt(exp1)

ggplot(exp1,aes(variable,value,fill=group))+
  geom_boxplot()+ 
  labs(x="",y="mRNA expression")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text=element_text(color='#333c41',size=12),
        legend.text = element_text(color='#333c41',size=12),
        legend.title = element_text(color='#333c41',size=13))+
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1.5))+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 1.5, alpha = 0.4,fill="white") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 1.5, alpha = 0.4,fill="#F1F3F3") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 1.5, alpha = 0.4,fill="#F1E2F6") +
  scale_fill_manual(values=c("#62d7f6","#F0b240"))+
  annotate('text', label = '*', x =1, y =1.35,  size =5,color="black")+
  annotate('text', label = '*', x =2, y =1.35, size =5,color="black")+
  annotate('text', label = '*', x =3, y =1.35, size =5,color="black")+
  ggtitle("GSE229750") 

##Immune infiltration
express_cpm
class(express_cpm)
express_cpm=as.data.frame(express_cpm)
colnames(id2symbol)
express_cpm$ENTREZID=rownames(express_cpm)
express_cpm=merge(express_cpm,id2symbol)
express_cpm=express_cpm[!duplicated(express_cpm$SYMBOL),]
rownames(express_cpm)=express_cpm$SYMBO
colnames(express_cpm)
express_cpm=express_cpm[,-c(1,14)]
express_cpm=express_cpm[,-c(2,4)]
library(IOBR)

#cibersort
im_cibersort<-deconvo_tme(eset=express_cpm,
                          method="cibersort",
                          arrays=F,
                          perm=1000
)
group1=group[-c(2,4),]
colnames(group1)
colnames(im_cibersort)[1]="geo_accession"
im_cibersort=merge(im_cibersort,group1)
colnames(im_cibersort)

library(ggplot2)
library(reshape)
library(tidyverse)

colnames(im_cibersort)[2:26]
colnames(im_cibersort)[2:26]=str_split(colnames(im_cibersort)[2:26],"_CIBERSORT",simplify = T)[,1]
colnames(im_cibersort)
im_cibersort1=im_cibersort[,-c(1,24,25,26,27)]
im_cibersort1$group=c(rep("FM",5),rep("Control",5))

im_cibersort1=im_cibersort1 %>%
  pivot_longer(cols =!c("group") , names_to = "names", values_to = "values")

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(), 
                 axis.text.x = element_text(angle = 90, hjust = 1 ),
                 panel.grid=element_blank(), 
                 legend.position = "top", 
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

p1 <- ggplot(im_cibersort1, aes(x = names, y = values,fill=group))+  
  labs(y="Cell composition",x= NULL,title = "GSE229750")+    
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+  
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+  
  theme_classic() + mytheme +   
  stat_compare_means(aes(group =  group), 
                     label = "p.signif",     
                     method = "wilcox.test",  
                     hide.ns = T)
p1

save(file="20240701GSE229750.Rdata")
