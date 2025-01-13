##Raw data processing
rm(list = ls())
library(GEOquery)
gse_number = "GSE67311"
eSet <- getGEO(gse_number, destdir = '.', getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
exp <- exprs(eSet)
dim(exp)
exp[1:4,1:4]

#Clinical information
pd <- pData(eSet)
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
gpl_number <- eSet@annotation;gpl_number
save(file="step1output.Rdata")

# Group
library(stringr)
Group=ifelse(str_detect(pd$`source_name_ch1`,"Control"),
               "control",
               "FM")
Group = factor(Group,levels = c("control","FM"))
Group
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL11532')

#PCA
dat=as.data.frame(t(exp))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Group, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
pca_plot


#DEG
library(limma)
design=model.matrix(~Group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
ids = ids[!duplicated(ids$symbol),]
deg <- inner_join(deg,ids,by="probe_id")
nrow(deg)

logFC_t=0
P.Value_t = 0.05
table(deg$change)
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
save(Group,deg,logFC_t,P.Value_t,gse_number,file = "step4output.Rdata")
write.csv(deg,"GSE67311deg.csv")

##Expressing bar chart
###7909225  DYRK3
#8130394 RGS17
#8109161 ARHGEF37

which(rownames(exp)==7909225) ##5437
which(rownames(exp)==8130394) ##27617
which(rownames(exp)==8109161) ##25420

exp1=exp[c(5437,27617,25420),]
rownames(exp1)=c("DYRK3","RGS17","ARHGEF37")
exp1=as.data.frame(t(exp1))
exp1$group=Group

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
  scale_y_continuous(expand = c(0, 0), limit = c(0, 10))+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 10, alpha = 0.4,fill="white") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 10, alpha = 0.4,fill="#F1F3F3") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 10, alpha = 0.4,fill="#F1E2F6") +
  scale_fill_manual(values=c("#62d7f6","#F0b240"))+
  annotate('text', label = '*', x =1, y =9,  size =5,color="black")+
  annotate('text', label = '****', x =2, y =9, size =5,color="black")+
  annotate('text', label = '*', x =3, y =9, size =5,color="black")+
  ggtitle("GSE67311") 

##Immune infiltration
library(IOBR)

class(exp)
exp=as.data.frame(exp)
colnames(ids)
exp$probe_id=rownames(exp)
exp=merge(exp,ids)
exp=exp[!duplicated(exp$symbol),]
rownames(exp)=exp$symbol
colnames(exp)
exp=exp[,-c(1,144)]

im_cibersort<-deconvo_tme(eset=exp,
                          method="cibersort",
                          arrays=F,
                          perm=1000
)
Group
im_cibersort$group=Group


library(ggplot2)
library(reshape)
library(tidyverse)

colnames(im_cibersort)[2:26]
colnames(im_cibersort)[2:26]=str_split(colnames(im_cibersort)[2:26],"_CIBERSORT",simplify = T)[,1]
colnames(im_cibersort)
im_cibersort1=im_cibersort[,-c(1,24,25,26)]

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
  labs(y="Cell composition",x= NULL,title = "GSE67311")+    
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+  
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+  
  theme_classic() + mytheme +   
  stat_compare_means(aes(group =  group), 
                     label = "p.signif",     
                     method = "wilcox.test",  
                     hide.ns = T)
p1

#Volcanic map
library(dplyr)
library(ggplot2)
dat  = deg[!duplicated(deg$symbol),]

p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

for_label <- dat%>% 
  filter(symbol %in% c("DYRK3","RGS17","ARHGEF37"))

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )

volcano_plot
