
rm(list = ls())
load("IMgene.Rdata")

############RGS17#######################
RGS17_IMGENE
library(BioEnricher)
all.enrich <- lzq_ORA.integrated(
  genes = RGS17_IMGENE,
  background.genes = NULL,
  GO.ont = 'ALL',
  perform.WikiPathways = T,
  perform.Reactome = T,
  perform.MsigDB = T,
  MsigDB.category = 'ALL',
  perform.Cancer.Gene.Network = T,
  perform.disease.ontoloty = T,
  perform.DisGeNET = T,
  perform.CellMarker = T,
  perform.CMAP = T,
  min.Geneset.Size = 3
)

lzq_ORA.barplot1(enrich.obj = all.enrich$simplyGO)
lzq_ORA.dotplot1(enrich.obj = all.enrich$simplyGO)

lzq_ORA.barplot1(enrich.obj = all.enrich$KEGG)
lzq_ORA.dotplot1(enrich.obj = all.enrich$KEGG)

lzq_ORA.barplot1(enrich.obj = all.enrich$WikiPathways)
lzq_ORA.dotplot1(enrich.obj = all.enrich$WikiPathways)

lzq_ORA.barplot1(enrich.obj = all.enrich$Reactom)
lzq_ORA.dotplot1(enrich.obj = all.enrich$Reactom)

lzq_ORA.barplot1(enrich.obj = all.enrich$DiseaseOntology)
lzq_ORA.dotplot1(enrich.obj = all.enrich$DiseaseOntology)

lzq_ORA.barplot1(enrich.obj = all.enrich$DisGeNET)
lzq_ORA.dotplot1(enrich.obj = all.enrich$DisGeNET)


lzq_ORA.barplot1(enrich.obj = all.enrich$GO)
lzq_ORA.dotplot1(enrich.obj = all.enrich$GO)

lzq_ORA.barplot1(all.enrich$Module.KEGG)

library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(viridis)

allgg <- as.data.frame(all.enrich$GO)
top5 <- allgg %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)
df_top5 <-rbind(subset(top5, ONTOLOGY=="BP"),subset(top5, ONTOLOGY=="CC"),subset(top5, ONTOLOGY=="MF"))
df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c("BP","CC",'MF'))
df_top5$Description <- factor(df_top5$Description, levels = rev(df_top5$Description))

options(repr.plot.width=7, repr.plot.height=6.5)

# mycol3 <- c('#6BA5CE', '#F5AA5F')
mycol3 <- c( '#BE88DC',"#6BA5CE","#F5AA5F")
cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")

p <- ggplot(data = df_top5, aes(x = Count, y = Description, fill=ONTOLOGY)) +
  geom_bar(width = 0.5,stat = 'identity') +
  theme_classic() + 
  scale_x_continuous(expand = c(0,0.5)) +
  scale_fill_manual(values = alpha(mycol3, 0.66))

p <- p + theme(axis.text.y = element_blank()) + 
  geom_text(data = df_top5,
            aes(x = 0.1, y = Description, label = Description),
            size = 4.8,
            hjust = 0) 

p <- p + geom_text(data = df_top5,
                   aes(x = 0.1, y = Description, label = geneID , color=-log10(p.adjust)),
                   size = 3.5,
                   fontface = 'italic',
                   hjust = 0,
                   vjust = 3) +
  scale_colour_viridis(option=cmap[7], direction=-1) 

p <- p + labs(title = 'RGS17 Enriched top5 GO') +
  theme(
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank())

p
ggsave("RGS17_gene_GO.pdf", p, w=11, h=9)

########all.enrich$simplyGO
allgg <- as.data.frame(all.enrich$simplyGO)
top5 <- allgg %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5)
df_top5 <-rbind(subset(top5, ONTOLOGY=="BP"),subset(top5, ONTOLOGY=="CC"),subset(top5, ONTOLOGY=="MF"))
df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c("BP","CC",'MF'))
df_top5$Description <- factor(df_top5$Description, levels = rev(df_top5$Description))

options(repr.plot.width=7, repr.plot.height=6.5)

# mycol3 <- c('#6BA5CE', '#F5AA5F')
mycol3 <- c( '#BE88DC',"#6BA5CE","#F5AA5F")
cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")

p <- ggplot(data = df_top5, aes(x = Count, y = Description, fill=ONTOLOGY)) +
  geom_bar(width = 0.5,stat = 'identity') +
  theme_classic() + 
  scale_x_continuous(expand = c(0,0.5)) +
  scale_fill_manual(values = alpha(mycol3, 0.66))

p <- p + theme(axis.text.y = element_blank()) + 
  geom_text(data = df_top5,
            aes(x = 0.1, y = Description, label = Description),
            size = 4.8,
            hjust = 0) 

p <- p + geom_text(data = df_top5,
                   aes(x = 0.1, y = Description, label = geneID , color=-log10(p.adjust)),
                   size = 4,
                   fontface = 'italic',
                   hjust = 0,
                   vjust = 3) +
  scale_colour_viridis(option=cmap[7], direction=-1) 

p <- p + labs(title = 'RGS17 Enriched top 5 simplyGO') +
  theme(
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank())

p
ggsave("RGS17_gene_simplyGO.pdf", p, w=11, h=9)



