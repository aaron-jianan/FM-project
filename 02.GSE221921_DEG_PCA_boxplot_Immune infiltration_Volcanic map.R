rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(data.table)

#1.Raw data processing
library(readxl)
counts<- read_excel("GSE221921_FM_ProcessedData.xlsx", 
                    sheet = "Values (FPKM)")
clinical <-  read_excel("GSE221921_FM_ProcessedData.xlsx", 
                        sheet = "Metadata (Samples)")
colnames(counts)
rawcount=counts[,c(2,9:ncol(counts))]
class(rawcount)
str(rawcount)
rawcount=as.data.frame(rawcount)
rawcount=rawcount[!duplicated(rawcount$Hugo_Gene_Symbol),]
names1=rawcount$Hugo_Gene_Symbol
rownames(rawcount)=names1
rawcount=rawcount[,-1]
rawcount=log2(rawcount+1)
rawcount[1:4,1:4]

keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
table(keep)
filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)

# cpm
library(edgeR)
express_cpm <- cpm(filter_count)
express_cpm[1:6,1:6]

##PCA
library(FactoMineR)
library(factoextra)
library(corrplot)
library(pheatmap)
library(tidyverse)

group <- clinical
colnames(group)

group <- group[match(colnames(rawcount), group$Sample),]
group

dat <- log10(express_cpm+1)
dat[1:4,1:4]
dim(dat)
sampleTree <- hclust(dist(t(dat)), method = "average")
plot(sampleTree)

temp <- as.data.frame(cutree(sampleTree,k = 2)) %>% 
  rownames_to_column(var="sample")
temp1 <- merge(temp,group,by.x = "sample",by.y="Sample")
table(temp1$`cutree(sampleTree, k = 2)`,temp1$Etiology)

dat <- log10(express_cpm+1)
dat[1:4,1:4]
dat <- as.data.frame(t(dat))
dat_pca <- PCA(dat, graph = FALSE)
dat[1:4,1:4]

group_list <- group[match(group$Sample,rownames(dat)),2]
group_list$Etiology=as.factor(group_list$Etiology)

mythe <- theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

p <- fviz_pca_ind(dat_pca,
                  geom.ind = "text", 
                  col.ind = group_list$Etiology, 
                  palette = c("#00AFBB", "#E7B800"), 
                  addEllipses = T,  
                  legend.title = "Groups") + mythe
p

##DEG
library(edgeR)
library(ggplot2)

filter_count[1:4,1:4]
group <- clinical

group_list <- group[match(colnames(filter_count),group$Sample),2]
group_list$Etiology=as.factor(group_list$Etiology)
table(group_list)

design <- model.matrix(~0+group_list$Etiology)
rownames(design) <- colnames(filter_count)
colnames(design) <- levels(factor(group_list$Etiology))
design

DEG <- DGEList(counts=filter_count, 
               group=factor(group_list$Etiology))

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
DEG_edgeR$GeneID=rownames(DEG_edgeR)


library(tidyverse)
DEG_edgeR_Sig <- filter(DEG_edgeR,regulated!="normal")

write.csv(DEG_edgeR,"4.DEG_edgeR_all.csv", row.names = F)
write.csv(DEG_edgeR_Sig,"4.DEG_edgeR_Sig.csv", row.names = F)
save(DEG_edgeR,DEG_edgeR_Sig,file = "Step03-edgeR_nrDEG.Rdata")


##Box diagram representation
###7909225  DYRK3
#8130394 RGS17
#8109161 ARHGEF37

which(rownames(express_cpm)=="DYRK3") ##7593
which(rownames(express_cpm)=="RGS17") ##1825
which(rownames(express_cpm)=="ARHGEF37") ##13354

exp1=express_cpm[c(7593,1825,13354),]
rownames(exp1)=c("DYRK3","RGS17","ARHGEF37")
exp1=as.data.frame(t(exp1))
exp1$group=group$Etiology

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
  scale_y_continuous(expand = c(0, 0), limit = c(0, 30))+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 30, alpha = 0.4,fill="white") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 30, alpha = 0.4,fill="#F1F3F3") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 30, alpha = 0.4,fill="#F1E2F6") +
  scale_fill_manual(values=c("#62d7f6","#F0b240"))+
  annotate('text', label = '***', x =1, y =29,  size =5,color="black")+
  annotate('text', label = '***', x =2, y =29, size =5,color="black")+
  annotate('text', label = '****', x =3, y =29, size =5,color="black")+
  ggtitle("GSE221921") 

##Immune infiltration
library(IOBR)
express_cpm

##cibersort
im_cibersort<-deconvo_tme(eset=express_cpm,
                          method="cibersort",
                          arrays=F,
                          perm=1000
)
colnames(im_cibersort)[1]="Sample"
im_cibersort=merge(im_cibersort,group)
colnames(im_cibersort)

library(ggplot2)
library(reshape)
library(tidyverse)

colnames(im_cibersort)[2:26]
colnames(im_cibersort)[2:26]=str_split(colnames(im_cibersort)[2:26],"_CIBERSORT",simplify = T)[,1]
colnames(im_cibersort)
im_cibersort1=im_cibersort[,-c(1,24,25,26,28)]

im_cibersort1=im_cibersort1 %>%
  pivot_longer(cols =!c("Etiology") , names_to = "names", values_to = "values")

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

p1 <- ggplot(im_cibersort1, aes(x = names, y = values,fill=Etiology))+  
  labs(y="Cell composition",x= NULL,title = "GSE221921")+    
  geom_boxplot(aes(fill = Etiology),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+  
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+  
  theme_classic() + mytheme +   
  stat_compare_means(aes(group =  Etiology), 
                     label = "p.signif",     
                     method = "wilcox.test",  
                     hide.ns = T)
p1

##Volcanic map
library(ggplot2)
library(tidyverse)
data <- DEG_edgeR
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
which(rownames(data)=="DYRK3")#2879
which(rownames(data)=="RGS17")#3212
which(rownames(data)=="ARHGEF37")#439

label <- data[c(2879,439,3212),]

p1 <- p + geom_point(size = 1, shape = 2, data = label,color="black") +
  ggrepel::geom_text_repel( aes(label = GeneID), data = label, color="black" )

p1
png(file = "5.Volcano_Plot.png",width = 900, height = 800, res=150)
plot(p1)
dev.off()

save(file = "20240701GSE221921.Rdata)