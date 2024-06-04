#### Script4 #####
# Obj: RNA-seq of KIT and MH63 inoculated with or without NSM , including: 
# Differential analysis using DESeq2
  # Fig3a-b. Volcano plot 
  # Fig3c. heatmap
  # Fig3d. GO enrichment visualization
  # Fig3e. KEGG enrichment visualization


library(tximport)
library(pheatmap)
source("save_pheatmap.R")

## 1. Load TPM results from Salmon
path<- "~/Documents/Research/LabMember/xjw/Publish/DataAvailibility/Code/"
setwd(path)

tx2gene<- read.csv("./source/4_RNAseq_Kit_and_MH63_NSM/tx2gene.csv",header = T)
dim(tx2gene)
head(tx2gene)
path<- "source/4_RNAseq_Kit_and_MH63_NSM/"
quant.sf<- file.path("source/4_RNAseq_Kit_and_MH63_NSM","quant",list.files("source/4_RNAseq_Kit_and_MH63_NSM/quant"))
names(quant.sf)<- gsub(".quant.sf","",list.files("source/4_RNAseq_Kit_and_MH63_NSM/quant"))
txi <- tximport(quant.sf, type="salmon", 
                tx2gene=tx2gene[,c("Name", "GeneID")], 
                countsFromAbundance="lengthScaledTPM")
names(txi)
txi_abund<- txi$abundance
txi_count<- txi$counts

# correlation between them:
sampabd_cor<- round(cor(txi_abund),digits = 2)
sampct_cor<- round(cor(txi_count),digits=2)

## Heatmap
sampabd_pdf<- pheatmap(sampabd_cor, display_numbers = T,fontsize = 10, angle_col = 45)
save_pheatmap_pdf(sampabd_pdf,"Fig3_Sup_sample_abundance_cor.pdf")
sampct_pdf<- pheatmap(sampct_cor, display_numbers = T,fontsize = 10, angle_col = 45)
save_pheatmap_pdf(sampct_pdf,"Fig3_Sup_sample_counts_cor.pdf")

## scatter plot
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DESeq2)
# Design
design<- read.csv("./source/4_RNAseq_Kit_and_MH63_NSM/Table2.Design.csv",header = T,row.names = 1)
design$Treatment<- factor(design$Treatment,levels = c("germfree","microbiota"))
#dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)


## MH63 DESEQ
#Tximport
subquant.sf<- quant.sf[grepl("MH63",basename(file.path(path,"quant",list.files("source/4_RNAseq_Kit_and_MH63_NSM/quant"))))]
subtxi <- tximport(subquant.sf, type="salmon", 
                   tx2gene=tx2gene[,c("Name", "GeneID")], 
                   countsFromAbundance="lengthScaledTPM")
# Design
MHdesign<- design[design$Cultivar=="MH63",]
MHdds <- DESeqDataSetFromTximport(subtxi, MHdesign, ~Treatment)
MHdds

# filter
keep <- rowSums(counts(MHdds)>0)>=2 & rowSums(counts(MHdds)) > 20
table(keep)
MHdds <- MHdds[keep, ]
MHdds <- DESeq(MHdds)
summary(results(MHdds))

res<- results(MHdds)

res <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res <- res[order(res$pvalue, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res <- cbind(subtxi$abundance[rownames(res),],res)
write.csv(res,"Fig3_Sup_MH63_res_all.csv")

LFC=1
pvalue=0.05
res_up<- res[which(res$log2FoldChange >= LFC & res$padj < pvalue),]      
res_down<- res[which(res$log2FoldChange <= -1*LFC & res$padj < pvalue),]   
res_total <- rbind(res_up,res_down)
write.csv(res_total,"Fig3_Sup_MH63_sig_all.csv")



## Kitaake DESEQ
## Kit
#Tximport
subquant.sf<- quant.sf[grepl("kit",basename(file.path(path,"quant",list.files("source/4_RNAseq_Kit_and_MH63_NSM/quant"))))]
subtxi <- tximport(subquant.sf, type="salmon", 
                   tx2gene=tx2gene[,c("Name", "GeneID")], 
                   countsFromAbundance="lengthScaledTPM")
Kitdesign<- design[design$Cultivar=="Kitaake",]
Kitdds <- DESeqDataSetFromTximport(subtxi, Kitdesign, ~Treatment)
Kitdds
# filter
keep <- rowSums(counts(Kitdds)>0)>=3 & rowSums(counts(Kitdds)) > 20
table(keep)
Kitdds <- Kitdds[keep, ]
Kitdds <- DESeq(Kitdds)
summary(results(Kitdds))

res<- results(Kitdds)

res <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res <- res[order(res$pvalue, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res <- cbind(subtxi$abundance[rownames(res),],res)
write.csv(res,"Fig3_Sup_Kit_res_all.csv")

LFC=1
pvalue=0.05
res_up<- res[which(res$log2FoldChange >= LFC & res$padj < pvalue),]      
res_down<- res[which(res$log2FoldChange <= -1*LFC & res$padj < pvalue),]   
res_total <- rbind(res_up,res_down)

write.csv(res_total,"Fig3_Sup_Kit_sig_all.csv")

###### compared with edgeR
# load("")
# kitup<- is.DEG_Kit.MvsG[is.DEG_Kit.MvsG==1,] %>% rownames()
# kitdown<- is.DEG_Kit.MvsG[is.DEG_Kit.MvsG==-1,] %>% rownames()
# kitdeg<- is.DEG_Kit.MvsG[is.DEG_Kit.MvsG!=0,] %>% rownames()


## KIT volcano plot (Fig3a) ####
Kit_res<- read.csv("Fig3_Sup_Kit_res_all.csv",header = T,row.names = 1)
Kit_res

#Add change_label
Res<- Kit_res
Res$Change<- "NS"
Res$Change[Res$log2FoldChange> 1 & Res$padj<0.05]<- "Up"
Res$Change[Res$log2FoldChange< -1 & Res$padj<0.05]<- "Down"
Kit_res$Change<- Res$Change

Res<- Kit_res
table(Res$Change)
ggplot(data = Res, aes(x = log2FoldChange, y = -log10(padj), col = Change)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2,alpha=.5) +
  theme_bw()+
  theme(axis.text = element_text(size=15))+
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Down", "NS", "Up")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  #coord_cartesian(ylim = c(0, 60), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),
       caption="Down:2540; Up: 1608; NS:20277") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) + # to customise the breaks in the x axis
  ggtitle('Kitaake Root: Microbiota_vs_Germfree')
ggsave("Fig4a.Kitaake_volcano.pdf",width = 12.5,height = 10,units = "cm")



## MH63 volcano plot (Fig3a) ####
MH_res<- read.csv("Fig3_Sup_MH63_res_all.csv",header = T,row.names = 1)
MH_res

#Add change_label
Res<- MH_res
Res$Change<- "NS"
Res$Change[Res$log2FoldChange> 1 & Res$padj<0.05]<- "Up"
Res$Change[Res$log2FoldChange< -1 & Res$padj<0.05]<- "Down"
MH_res$Change<- Res$Change

Res<- MH_res
table(Res$Change)
ggplot(data = Res, aes(x = log2FoldChange, y = -log10(padj), col = Change)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2,alpha=.5) +
  theme_bw()+
  theme(axis.text = element_text(size=15))+
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Down", "NS", "Up")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  #coord_cartesian(ylim = c(0, 60), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Severe', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),
       caption="Down:1987; Up: 2059; NS:18187") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) + # to customise the breaks in the x axis
  ggtitle('MH63 Root: Microbiota_vs_Germfree')
ggsave("Fig4b.MH63_volcano.pdf",width = 12.5,height = 10,units = "cm")

## GO plot (Fig3d) ####
## GO
# The readin file(Table2.GO_BP_sig0.01.txt) was modified from Version1(`"Version1/4.GO/MH_sig_GO.xls"` and `"Version1/4.GO/Kit_sig_GO.xls"`). Only BP class with padj <0.01 was retained.
# path<- "/Users/hupeiyao/Documents/Research/LabMember/xjw/RNA-seq/20231227/Version4/";setwd(path)
gobp<- read.delim("source/4_RNAseq_Kit_and_MH63_NSM/Table2.GO_BP_sig0.01.txt",sep="\t",header = T)

gobp$Type<- factor(gobp$Type)
gobp$xlab<- paste(gobp$Type,gobp$GO_Name,sep=":")
gobp$padj_log10<- -log10(gobp$padj+10^-20)
gobp$Term<- paste0(gobp$GO_ID,": ",gobp$GO_Name)
golevel<- c("GO:0006810: transport" ,"GO:0009607: response to biotic stimulus","GO:0009719: response to endogenous stimulus",
            "GO:0009628: response to abiotic stimulus","GO:0006950: response to stress"  ,"GO:0050896: response to stimulus",
            "GO:0019748: secondary metabolic process","GO:0008152: metabolic process")


gobp$Term<- factor(gobp$Term,levels = golevel)

gobp %>% 
  ggplot(aes(x=Term,y=GeneNo,fill=padj_log10,color=Type))+ 
  geom_bar(stat = "identity", width = .5,position=position_dodge(width=.65),size=.3)+    theme_bw()+
  coord_flip()+
  labs(x="Term",y = "RichFactor",title = "DEGs GO-BP")+ 
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(color="black"))+
  scale_fill_gradient(low = "#F2E8D4", high = "#8E9F85")+
  #scale_fill_gradient(low = "#D3D3D3", high = "#7A7A7A")+
  scale_color_manual(values=c("#E1522F","#1E7CAE"))+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())
ggsave("Fig3d.GO_col2.pdf",height = 3,width = 7.5)


## KEGG

kegg_all<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/KEGG_merge.txt",header = T,sep="\t")
names(kegg_all)[3]<- "GeneNo"
names(kegg_all)[10]<- "padj"

data1<- kegg_all %>%
  dplyr::select(Term.Name,GeneNo) %>%
  group_by(Term.Name) %>%
  summarise_each(funs=mean) %>%
  dplyr::arrange(GeneNo)
levelkegg<- data1$Term.Name
kegg_all$Term.Name[order(kegg_all$GeneNo)]

#kolevel<- dplyr::arrange(kegg_all,desc(GeneNo))$Term.Name %>% unique
#kolevel<- dplyr::arrange(kegg_all,padj)$Term.Name %>% unique
kegg_all$Term.Name<- factor(kegg_all$Term.Name,levels = levelkegg)
kegg_all$padj_log10<- -log10(kegg_all$padj+10^-20)

kegg_all %>%
  ggplot(aes(x=Term.Name,y=GeneNo,fill=padj_log10,color=Type))+ 
  geom_bar(stat = "identity", width = .5,position=position_dodge(width=.6),size=.3)+    theme_bw()+
  coord_flip()+
  labs(x="Term",y = "RichFactor",title = "GO0008152 KEGG Enrichment")+ 
  scale_fill_gradient(low = "#F2E8D4", high = "#8E9F85")+
  #scale_fill_gradient(low = "#D3D3D3", high = "#7A7A7A")+
  scale_color_manual(values=c("#E1522F","#1E7CAE"))+  
  theme(#axis.text.y=element_text(angle=30),
    #axis.text.y = element_text(size=20),
    axis.text = element_text(color = "black"),
    axis.title = element_text(colour = "black"))

ggsave("Fig3e.KEGG.pdf",width = 8,height = 5)
write.csv(kegg_all,"Fig3e.KEGG_plot.csv")


### Figure3c #####

library('readr')
library('magrittr')
library('tibble')
library('gplots')
library('dendextend')
library('dynamicTreeCut')
library('ggplot2')
library('tidyr')
library('DESeq2')
library('dplyr')
library('RColorBrewer')
library('gridExtra')
library('cluster')
library('scales')


######## Construction #########
dds <- DESeqDataSetFromTximport(txi, design, ~Treatment)
keep <- rowSums(counts(dds)>0)>=5 & rowSums(counts(dds)) > 20
table(keep)
dds <- dds[keep, ]
dds <- DESeq(dds)
summary(results(dds))

####### View normalization #########
colors<- c(rep("#023C88",3),rep("#52B2CF",3),
           rep("#527C5A",2),rep("#8E9F85",2))
par(mfrow = c(2,1))
pdf("Fig3_Sup_normalize.pdf",width=10,height = 5)
boxplot(log10(counts(dds,normalized=F)),main="Unnormalized abundance", ylab="log10(abundance)",col=colors)
boxplot(log10(counts(dds,normalized=T)),main="Normalized abundance", ylab="log10(abundance)",col=colors)
dev.off()

## count transformation
# count transformation #######
rld <- rlog(dds)   #class: DESeqTransform 
ntd <- normTransform(dds)
dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

## remove genes with zero expression
dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}
ddssvaRmZero <- ddssva[rownames(ddssva) %in% rownames(dat), ]

ddssvaRmZero <- DESeq(ddssvaRmZero)

resultsNames(ddssvaRmZero)

expre<- data.frame(counts(ddssvaRmZero,normalized=T))
results_treat<- results(ddssvaRmZero) %>% data.frame
results_treat<- cbind(expre[rownames(results_treat),],results_treat)



## Data setup
library(dplyr)
library(magrittr)
library(tibble)
# significant genes
Kit_sig<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/Table1.Kit_sig_all_norm.csv",header = T,row.names = 1)
Kit_sig<- Kit_sig[grep("kit",colnames(Kit_sig))]
MH_sig<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/Table1.MH63_sig_all_norm.csv",header = T,row.names = 1)
MH_sig<-MH_sig[grep("MH",colnames(MH_sig))]

kit_res<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/Table1.Kit_res_all_norm.csv",header = T,row.names = 1)
MH_res<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/Table1.MH63_res_all_norm.csv",header = T,row.names = 1)
kit_res$LOCUS<- rownames(kit_res)
MH_res$LOCUS<- rownames(MH_res)


sig_loc<- c(rownames(Kit_sig),rownames(MH_sig)) %>% unique()
length(sig_loc)  #6643

## stats info:
idx<- c("log2FoldChange","pvalue","padj","LOCUS")
res_stat<- cbind(kit_res[sig_loc,idx],MH_res[sig_loc,idx]);dim(res_stat)
head(res_stat)
res_stat<- res_stat[,-grep("LOCUS",colnames(res_stat))]
colnames(res_stat)=c(paste0("Kit_microbiota_vs_germfree_",idx[1:3]),
                     paste0("MH63_microbiota_vs_germfree_",idx[1:3]))
res_stat %<>% replace(is.na (.), 0)
rownames(res_stat)<- sig_loc
###

# merge Kit and MH, and substitute NA with 0 (low expression genes removed in one of the cultivar)
all_sig<- cbind(kit_res[sig_loc,],MH_res[sig_loc,])
sampleN=colnames(all_sig)[c(1:6,14:17)]

all_sig_to0<- all_sig %>% replace(is.na (.), 0)
all_sig_to0<- all_sig_to0[,sampleN]

head(all_sig_to0)
head(all_sig)

#
ncbi_conversion<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/NCBI_conversion_final.txt",sep="\t",header = T)
all_sig_info<- cbind(all_sig_to0,res_stat[rownames(all_sig_to0),])
write.csv(all_sig_info,"Fig3c_Sup.all_sig_to0_count.csv")


## k-means clustering
z_var <- apply(all_sig_to0, 1, var)
z_mean <- apply(all_sig_to0, 1, mean)
pdf("Fig3c_Sup.z_mean_var.pdf",width = 5,height = 5)
plot(log2(z_mean), log2(z_var), pch = '.',xlim=c(0,25),ylim=c(0,35))
abline(h = 1, col='red')
abline(v = 1, col='red')
text(x = 5,
     y = 23,
     labels = 'variance > 1 &\n mean > 1',
     col = 'red')
dev.off()

#scaling
scaleCount<- all_sig_to0 %>% t() %>% scale() %>% t()
summary(scaleCount[1,])
dim(scaleCount)
head(scaleCount)
table(scaleCount %>% complete.cases())

########### Tips ##############
# check for complete cases, i.e. rows without missing values. 
# It returns  a logical vector  indicating  whether each row is complete or not.
library(magrittr) #`%<>%`: This is the magrittr pipe assignment operator
a<- all_sig
a %<>% .[complete.cases(.),]
table(complete.cases(all_sig))
dim(a)
#############################

######
# `Elbow Method` to determine optimal clusters
# ref: https://uc-r.github.io/kmeans_clustering
library(factoextra) # clustering algorithms & visualization
set.seed(123)
pdf("Fig3c_Sup_Figure4.kmeans_sse.pdf",width = 5,height = 5)
fviz_nbclust(scaleCount,kmeans,method = "wss")
dev.off()
#######

## execute
set.seed(26)
kClust10 <- kmeans(scaleCount, centers = 10, algorithm= 'MacQueen', nstart = 1000, iter.max = 20)
kClust6 <- kmeans(scaleCount, centers = 6, algorithm= 'MacQueen', nstart = 1000, iter.max = 20)
kClust5 <- kmeans(scaleCount, centers = 5, algorithm= 'MacQueen', nstart = 1000, iter.max = 20)
kClust4 <- kmeans(scaleCount, centers = 4, algorithm= 'MacQueen', nstart = 1000, iter.max = 20)

table(kClust10$cluster)
table(kClust6$cluster)
table(kClust4$cluster)

## selected genes
cl <- kClust10$cluster
prefix <- 'kmeans_10'

clusterGene <- scaleCount %>%
  as.data.frame %>%
  rownames_to_column(var = 'ID') %>%
  as_tibble %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  }

## plot core cluster
clusterCore <- clusterGene %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', cl))) %>%
  gather(Sample, NorExpress, -1) %>%
  mutate(Sample = Sample )
# write.csv(clusterCore,"Fig3d_Sup.kmeans_corecluster.csv")

## plot all genes
clusterGenePlot <- clusterGene %>%
  gather(Sample, NorExpress, -ID, -cl) %>%
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', sort(unique(cl))))) %>%
  mutate(Sample = Sample)

clusterCorePlot <- clusterCore %>% dplyr::mutate(ID = 1 : nrow(clusterCore))
ggplot(clusterGenePlot, aes(Sample, NorExpress, group = ID)) +
  geom_line(color = 'grey30', alpha = 0.01) +
  facet_wrap(. ~ cl, ncol = 2) +
  geom_point(data = clusterCorePlot, aes(Sample, NorExpress, col = cl, group = ID)) +
  geom_line(data = clusterCorePlot, aes(Sample, NorExpress, group = cl, col = cl)) +
  ylab('Scaled counts') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(colour = guide_legend(title = 'kmeans (k=10)'))
ggsave(paste0("Fig3d_Sup_",prefix, '_genes.pdf'), width = 10, dpi = 320)
write.csv(clusterCorePlot,"Fig3d_Sup.clustercoreplot.csv")

## cluster cor pheno
##p value calculation from WGCNA
corPvalueStudent <- function(cor, nSamples) {
  
  ## ref: https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
  T <- sqrt(nSamples - 2) * cor / sqrt(1 - cor^2)
  
  p <- apply(T, 1:2, function(x) {
    if (x < 0) {
      eachp <- 1 -  pt(x, nSamples - 2, lower.tail = FALSE)
    } else {
      eachp <- pt(x, nSamples - 2, lower.tail = FALSE)
    }
  })
  
  return(p)
}
##~~~~~~~~~~~~~~~~~~~~~~~~cluster cor phenotype~~~~~~~~~~~~~~~~~
traits <- data.frame(Cultivar = c(0,0,0,0,0,0,1,1,1,1), #MH:1
                     Live = c(0,0,0,1,1,1,0,0,1,1)) #Microbiota: 1

cores <- clusterGene %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>%
  mutate(cl = cl %>% paste0('cluster_', .)) %>%
  column_to_rownames(var = 'cl') %>%
  t

moduleTraitCor <- cor(cores, traits, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(traits))
# write.csv(moduleTraitCor,"Fig3d_Sup.moduleTraitCor.csv")
# write.csv(moduleTraitPvalue,"Fig3d_Sup.moduleTraitPvalue.csv")

traitPPlot <- moduleTraitPvalue %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, pvalue, -1) %>%
  as_tibble

traitCorPlot <- moduleTraitCor %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, correlation, -1) %>%
  as_tibble %>%
  mutate(x = rep(0 : (ncol(traits) - 1), each = ncol(cores))) %>%
  mutate(y = rep((ncol(cores) - 1) : 0, ncol(traits))) %>%
  inner_join(traitPPlot) %>%
  mutate(addtext = paste0(round(correlation, digit = 2),
                          '\n',
                          '(',
                          round(pvalue, digit = 2),
                          ')'))

ggplot(traitCorPlot, aes(x = x, y = y, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100),
                       breaks = seq(-1, 1, 0.5),
                       labels = format(seq(-1, 1, 0.5)),
                       limits = c(-1, 1)) +
  geom_text(aes(label = addtext)) +
  scale_x_continuous(breaks = 0 : (ncol(traits) - 1), labels = colnames(traits)) +
  scale_y_continuous(breaks = 0 : (ncol(cores) - 1), labels = paste0('cluster_', (ncol(cores)):1)) +
  xlab('Trait') +
  ylab('Cluster')+theme_light()
ggsave(paste0("Fig3d_Sup.",prefix, '_trait.pdf'))


## Heatmap
scaleC <- all_sig_to0 %>%
  t %>%
  scale %>%
  t %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Scale_', .)))

rawC <- all_sig_to0 %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Raw_', .)))

degresC<- res_stat
degresC$ID<- rownames(degresC)

heatPlot <- rawC %>%
  inner_join(scaleC) %>%
  inner_join(degresC) %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  } %T>%
  {(sum(names(cl) == .$ID) == nrow(.)) %>% print} %>% ## check cl names and degresC row names
  dplyr::slice(cl %>% order)

deganno<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/Table4.res_stat2.csv",header = T,row.names = 1)
deganno$Gene<- rownames(deganno)
# cbind(heatPlot,deganno[heatPlot$ID,]) %>%
#   write_csv(paste0("../../Version2/Table/Table5.",prefix, '.csv'))

heatRawPlot <- heatPlot %>%
  dplyr::select(ID, starts_with('Raw')) %>%
  gather(sample, raw, -1) %>%
  mutate(x = rep(0 : 9, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 10))

heatScalePlot <- heatPlot %>%
  dplyr::select(ID, starts_with('Scale')) %>%
  gather(sample, scale, -1) %>%
  mutate(x = rep(0 : 9, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 10))

heatlog2FCPlot <- heatPlot %>%
  dplyr::select(ID, ends_with('FoldChange')) %>%
  gather(sample, log2FC, MH63_microbiota_vs_germfree_log2FoldChange : Kit_microbiota_vs_germfree_log2FoldChange) %>%
  mutate(x = rep(0 : 1, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 2))

## sig |FC| > 1 and padj < 0.05
fcsig <- heatPlot %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > 1 ~ 1,
                                 . < -1 ~ -1,
                                 TRUE ~ 0)))

padjsig <- heatPlot %>%
  dplyr::select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsigPlot <- (padjsig * fcsig) %>%
  as_tibble %>%
  gather(sample, sig, 1:2) %>%
  mutate(x = rep(0 : 1, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 2))

heatGroupPlot <- heatPlot %>%
  dplyr::select(ID, cluster = cl) %>%
  mutate(x = 0) %>%
  mutate(y = 0 : (nrow(heatPlot) - 1))

theme_bact <- function(...) {
  theme_bw() %+replace%
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks.length = unit(0, 'mm'),
          axis.line = element_blank(),
          panel.spacing = unit(0, 'mm'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), 'line'),
          legend.spacing = unit(0, 'mm'),
          ...)
}
sampleN=colnames(all_sig)[c(1:6,14:17)]

ggplot(heatScalePlot, aes(x = x, y = y, fill = scale)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral')))(100), name = 'scale(count)') +
  scale_x_continuous(breaks = 0 : 9,
                     labels = sampleN ) +
  theme_bact(legend.position = 'left',
             axis.text.x = element_text(angle = 90, hjust = 1))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~merge all plots~~~~~~~~~~~~~~~~~~~~
groupe <- ggplot(heatGroupPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(cluster))) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')

groupne <- heatGroupPlot %>%
  group_by(cluster) %>%
  summarise(y = median(y)) %>%
  mutate(x = 0, cluster = cluster %>% paste0('cluster', .)) %>%
  ggplot(aes(x = x, y = y, label = cluster)) +
  geom_text() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0),  limits = c(0, nrow(heatGroupPlot)), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')

rawe <- ggplot(heatRawPlot, aes(x = x, y = y, fill = log2(raw))) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name = 'GnBu'))(100), name = 'log2(count)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')

scalee <- ggplot(heatScalePlot, aes(x = x, y = y, fill = scale)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral')))(10), name = 'scale(count)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')

fce <- ggplot(heatlog2FCPlot, aes(x = x, y = y, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = 'PiYG')))(100), name = 'log2(FoldChange)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')

sige <- ggplot(heatsigPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(sig))) +
  scale_fill_manual(name = 'Significant', labels = c('Down', 'No', 'Up'), values = c('#023C88', '#A2A2A2', '#8F0D02')) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')

library(dplyr)
sigte <- heatGroupPlot %>%
  select(cluster, y) %>%
  inner_join(heatsigPlot) %>%
  select(sample, sig, x, y, cluster) %>%
  group_by(sample, cluster) %>%
  dplyr::count(sig) %>%
  spread(sig, n) %>%
  {
    loc <- heatGroupPlot %>%
      group_by(cluster) %>%
      summarise(y = median(y))
    inner_join(., loc)
  } %>%
  dplyr::rename('down' = `-1`, 'no' = `0`, 'up' = `1`) %>%
  ungroup %>%
  mutate_at(c('down', 'no', 'up'), .funs = list(~if_else(is.na(.), 0L, .))) %>%
  mutate(x = rep(c(0.2, 0.4, 0), each = max(cl))) %>%
  mutate(signum = paste0(down, '/', no, '/', up)) %>%
  select(signum, x, y) %>%
  ggplot(aes(x = x, y = y, label = signum)) +
  geom_text() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0),  limits = c(0, nrow(heatGroupPlot)), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.1, 0.5), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')



blanke <- ggplot(tibble(x = 0, y = 0 : (nrow(heatPlot) - 1)),
                 aes(x = x, y = y)) +
  geom_tile(colour = 'white') +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_bact(title = element_blank(),
             legend.position = 'none')

cairo_pdf(paste0("Fig3d_Sup.",prefix, '_heatmap_merge.pdf'), width = 15)
grid.arrange(groupne,
             groupe,
             blanke,
             rawe,
             blanke,
             scalee,
             blanke,
             fce,
             blanke,
             sige,
             
             nrow = 1,
             ncol = 11,
             widths = c(3.5, 1, 0.5, 13, 0.5, 13, 0.5, 3, 0.5, 3, 10) %>% {. / sum(.)})
dev.off()

g <- grid.arrange(groupne,
                  groupe,
                  blanke,
                  scalee,
                  nrow = 1,
                  ncol = 4,
                  widths = c(3.5, 1, 0.5, 13) %>% {. / sum(.)})
ggsave(file = paste0("Fig3d_Sup.",prefix, '_heatmap_all.pdf'), plot = g)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## write the cluster file
#inner_join(deganno, heatPlot) %>%
#  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
#  write_csv(paste0("../../Version2/Table/Table6.",prefix, '.csv'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################
```


## select DEGs


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('Fig3c.Sup_kmeans_10.csv')
kmeansRes <- wholeDEG %>%
  dplyr::select(ID, cl)

fcsig <- wholeDEG %>%
  dplyr::select(contains('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(2) ~ 1,
                                 . < -log2(2) ~ 1,
                                 TRUE ~ 0)))
padjsig <- wholeDEG %>%
  dplyr::select(contains('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

heatsig %>%
  mutate_at(c('ID', 'Product'), .funs = list(~if_else(is.na(.), '', .))) 
#%>%  write_csv('../../Version2/Table/Table5.heat_deg_sig.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


```

## Advanced plot
#Scaled counts normalized to library size were generated using DESeq2 (‘rlog’ function) and transformed as median-centered z- score (by transcripts, ‘scale’ function). Then, z-scores were used to conduct k-means clustering for all transcripts. The cluster number (k=10) was determined by sum of squared error and Akaike information criterion. Transcripts with similar expression pat- terns were grouped in the same cluster. Differentially expressed transcripts and cluster results were visualized using heatmaps generated with the ComplexHeatmap package in R (Gu et al., 2016). Gene expression in individual plots was plotted using scaled counts data. 

mhrld<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/Table3.MHdds_1_rlog.csv",row.names = 1)
kitrld<- read.csv("source/4_RNAseq_Kit_and_MH63_NSM/Table3.Kitdds_1_rlog.csv",row.names = 1)


## rlog transformed
rld<- cbind(mhrld[sig_loc,],kitrld[sig_loc,])
rownames(rld)<- sig_loc %>% replace(is.na(.),0)
rawC <- cbind(rld[heatsig$ID,],heatsig[,c("ID","cl")])

scaleC <- rawC[,1:10]%>%
  t %>% 
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% dplyr::select(ID, cl))
summary(as.numeric(scaleC[1,1:10]))


library(ComplexHeatmap)
cairo_pdf('Fig3c_kmeans_10_heatmap_sig2.pdf')
cultivar <- HeatmapAnnotation(Cultivar = c(rep("MH63",4),rep("Kitaake",6)),
                              col = list(Iron = c('Kitaake' = 'white', 'MH63' = 'grey')),
                              gp = gpar(col = 'black'))
treatment <- HeatmapAnnotation(Treatment = c(rep("Germfree",2),rep("Microbiota",2),
                                             rep("Germfree",3),rep("Microbiota",3)),
                               col = list(Treatment = c('Germfree' = 'white', 'Microbiota' = 'grey')),
                               gp = gpar(col = 'black'))
Heatmap(matrix = scaleC[,1:10],
        name = 'Scaled Counts',
        row_order = order(scaleC$cl) %>% rev,
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1:10,
        column_split = c(rep("MH63",4),(rep("Kitaake",6))),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -7, -8)])(10),
        top_annotation = c(cultivar, treatment))
dev.off()

write.csv(scaleC,"Fig3c.kmeans_10_heatmap_sig2_neo.csv")
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~box plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleN2 <- c('kit_germfree', 'kit_microbiota',
              'MH63_germfree', 'MH63_microbiota')
boxplotData<- rbind(apply(scaleC[,grep("kit_germfree",colnames(scaleC))],1,mean),
                    apply(scaleC[,grep("kit_microbiota",colnames(scaleC))],1,mean),
                    apply(scaleC[,grep("MH63_germfree",colnames(scaleC))],1,mean),
                    apply(scaleC[,grep("MH63_microbiota",colnames(scaleC))],1,mean)) %>% 
  t %>% 
  set_colnames(sampleN2) %>%
  set_rownames(scaleC$ID) 
boxplotData<- cbind(boxplotData[kmeansRes$ID,],kmeansRes)
head(boxplotData)




colors<- c(rep("#023C88",3),rep("#52B2CF",3),
           rep("#527C5A",2),rep("#8E9F85",2))
color2<- c("#023C88","#52B2CF","#527C5A","#8E9F85")
for (i in 1:10) {
  data1<- boxplotData %>%
    dplyr::filter(cl == i) %>%
    dplyr::select(-ID, -cl) %>%
    gather(key = 'Conditions', value = 'ScaleCounts') %>%
    mutate(Conditions = Conditions %>% factor(levels = sampleN2)) 
  data2<- data1 %>%
    group_by(Conditions) %>%
    summarise_each(funs=mean) %>% ungroup()
  
  data1 %>%
    ggplot(aes(x = Conditions, y = ScaleCounts,fill=Conditions),outlier.fill=Conditions) +
    geom_boxplot(lwd=.25,outlier.size=.25,width=.75) +
    geom_line(data=data2,mapping=aes(x=Conditions,y=ScaleCounts,group=1), color="#C70039")+ #linetype = "dashed",
    ylab(paste("Cluster",i,'Scaled counts')) +xlab("")+
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=18),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
          legend.text.align = 0,
          axis.line = element_line(size=.25),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.text=element_text(size= 13),
          legend.title = element_text(size = 14))+
    scale_color_manual("black")+scale_fill_manual(values = color2)
  #ggsave(paste0('../../Version2/Figure/Figure7/Figure7.kmeans_10_boxplot_cluster', i, '.pdf'),width = 15,height = 10,units = "cm")
  ggsave(paste0("Fig3c_kmeans_10_boxplot_cluster", i, '.pdf'),width = 15,height =7.5,units = "cm")
  
}

save.image("4_RNAseq_Kit_and_MH63_NSM.RData")
