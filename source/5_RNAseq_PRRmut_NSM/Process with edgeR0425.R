###############################################################################
## Xu Jiwei's RNA-seq data
##
##
##  2024.2.27

library(edgeR)
library("DESeq2")
library(ggplot2)
library(Glimma)
library(tximport)
library("cowplot")
library("ggvenn")
library("ComplexHeatmap")
library("RColorBrewer")
library(viridis);
library(ggfortify) #PCA plot
library(corrplot)
library("AnnotationHub")
library(UpSetR)

setwd("~/project/labmemeber/xjw/RNA-seq-XJW/RNA-seq-XJW2023")
OUT_DIR="./out_edgeR2/";

## Input the gene-level counts from salmon 
samples<-read.table("samples.txt", head=TRUE);
rownames(samples)<-samples$Seq_ID;
samples$condition<-as.factor(samples$condition);
samples$genotype<-as.factor(samples$genotype);

samples$group<-paste(samples$genotype, samples$condition, sep="_")
samples$group<- factor(samples$group, levels=c("kit_germfree", "kit_microbiota",
                                               "fls2_germfree", "fls2_microbiota",
                                               "cerk1_germfree", "cerk1_microbiota") )
samples<-samples[order(samples$group, decreasing = TRUE), ]
design <- model.matrix(~0+group, data = samples)
colnames(design) <- levels(samples$group)

files<-file.path("~/project/labmemeber/xjw/RNA-seq-XJW/RNA-seq-XJW2023/quants_salmon_NCBIref", 
                 paste(samples$Seq_ID, "quant", sep = "_"), "quant.sf");
length(files)
names(files)<-samples$Seq_ID;
file.exists(files);

tx2gene <- read.csv("~/project/pubdata/NIP_ncbi_dataset/data/GCF_001433935.1/NIP.tu2genes.csv", head=F);
txi <- tximport(files, type="salmon", tx2gene=tx2gene);
# the txi contains abundance, counts, length

# refRNA annotation was download from NCBI entrez FTP site, oryza sative
Os_entrez <- read.table("~/project/pubdata/NIP_ncbi_dataset/data/GCF_001433935.1/Oryza_sativa.gene_info", head=T, sep="\t",quote="");
Os_rapdb<- read.table("~/project/pubdata/NIP_ncbi_dataset/IRGSP-1.0_representative_annotation_2024-01-11.tsv", head=T, sep="\t",quote="",fill = TRUE)
Os_kegg<- read.table("~/project/pubdata/NIP_ncbi_dataset/kegg_list_osa.txt",head=T, sep="\t",quote="");
Os_funrice<-read.table("~/project/pubdata/NIP_ncbi_dataset/funrice_geneInfo.table.txt", header = T, sep="\t")
dim(Os_entrez) #35851    17
dim(Os_rapdb); #45019    19
length(unique(Os_rapdb$Locus_ID))  #37840
dim(Os_kegg)  # 28143
length(grep("OSNPB",Os_entrez$LocusTag))  #26123 has OSNPB ids
length(grep("-",Os_entrez$LocusTag))      #8208
table(Os_entrez$type_of_gene[(grep("-",Os_entrez$LocusTag))])

## Convert eh OSNPB tag to RAPDB tag
Os_entrez$RAPdb_ID<-"-";
t<-Os_entrez$LocusTag[grep("OSNPB",Os_entrez$LocusTag)]
Os_entrez$RAPdb_ID[grep("OSNPB",Os_entrez$LocusTag)]<-paste("Os", substr(t,7,8),"g",substr(t,9,15), sep="")
# Keep one RAPdb ID for each record. Note, some gene model have multiple RAPdb ID
Os_funrice$RAPdb<-substr(Os_funrice$RAPdb,1,15)
Os_funrice<-Os_funrice[grep("Os", Os_funrice$RAPdb),]
# Manually edit following records: ‘Os03g0793500’, ‘Os03g0819600’, ‘Os05g0200400’, ‘Os10g0524500’, ‘Os11g0547000’ 
rownames(Os_funrice)<-Os_funrice$RAPdb

Os_entrez$funrice<-"-"
t<-(Os_entrez$RAPdb_ID %in% Os_funrice$RAPdb)  #4140 characterized genes
Os_entrez$funrice[t]<-Os_funrice[Os_entrez$RAPdb_ID[t],1]

# calculate normalize factor following instructions from tximport
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for edgeR.
is_ref_rna<-stringr::str_detect(Os_entrez$Symbol, "LOC*")
Os_entrez<-Os_entrez[is_ref_rna,]
rownames(Os_entrez)<-Os_entrez$Symbol;

# Creating a DGEList object for edgeR.
y <- DGEList(cts, group = samples$group, genes = Os_entrez[rownames(cts),c(2:5,9:10,14,17:18)])
y$samples # check a assignment of group
y <- scaleOffset(y, normMat)

# filtering using the design information
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, ]
boxplot(log2(y$counts+1)) # without normalization
cpms <- edgeR::cpm(y, offset = y$offset, log = TRUE)
colnames(cpms)<-
  paste(samples[colnames(cpms),7], c("rep1","rep2","rep3"), sep="_")
boxplot(cpms)

# clustering samples. 
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255);
colors<-viridis(255)
sampleDists <- dist(t(cpms)); # the default method for distance measurement is "euclidean"  
sampleDistMatrix <- as.matrix( sampleDists )

p<-pheatmap(sampleDistMatrix,
            clustering_distance_rows=sampleDists,
            clustering_distance_cols=sampleDists,
            col=colors,show_colnames = TRUE,
            show_rownames = TRUE);

corr_samples<-cor(cpms)
rownames(corr_samples)
round(corr_samples,2)
pdf(file.path(OUT_DIR,"corrplot.pdf"))
corrplot(corr_samples, type="full", col.lim=c(0.9,1.0), is.corr=FALSE,
         order="hclust")
dev.off()

##MDS plot using plotMDS
p<-glimmaMDS(y,
          groups=data.frame(t(sapply(strsplit(as.character(y$samples$group), "_"),c))),
          top=30000)
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "RNA_seq_mdsplot.html"))

####
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
summary(fit$df.prior)

##########################################################################
## Check DEGs (microbiota vs germfree) in different genotypes
## microbiota vs germfree in Kit (WT), osfls2, and cerk1
## Comparable number of DEGs between EdgeR and DEseq2


Kit.MvsG <- makeContrasts(kit_microbiota-kit_germfree, levels=design)
res_Kit <- glmQLFTest(fit, contrast=Kit.MvsG)
is.DEG <- decideTestsDGE(res_Kit, p.value=0.05, lfc=1)
summary(is.DEG)

glimmaMA((res_Kit))
glimmaVolcano(res_Kit)
topTags(res_Kit)
res_Kit$table$logFC


fls2.MvsG <- makeContrasts(fls2_microbiota-fls2_germfree, levels=design)
res_microbe_fls2 <- glmQLFTest(fit, contrast=fls2.MvsG )
#res_microbe_fls2 <- glmTreat(fit, contrast = fls2.MvsG, lfc=log2(2))

glimmaMA(res_microbe_fls2)
glimmaVolcano(res_microbe_fls2)
topTags(res_microbe_fls2)
is.DEG <- decideTestsDGE(res_microbe_fls2, lfc = log2(2))
summary(is.DEG)


cerk.MvsG <- makeContrasts(cerk1_microbiota-cerk1_germfree, levels=design)
res_microbe_cerk <- glmQLFTest(fit, contrast = cerk.MvsG)
is.DEG <- decideTestsDGE(res_microbe_cerk, lfc = log2(2))
summary(is.DEG)

## The above DEG analysis show substantial variations of DEG numbers in WT, fls2, cerk1

########################################################################################
##  Genotype orient DEGs (mut vs WT) in two conditions

## DEGs in fls2 mutant vs WT in two conditions


fls2vsKit_germfree <- makeContrasts(fls2_germfree-kit_germfree, levels=design)
res_fls2vsKit_germfree <- glmQLFTest(fit, contrast=fls2vsKit_germfree )
hist(res_fls2vsKit_germfree$table$PValue)
topTags(res_fls2vsKit_germfree, sort.by = "logFC");
is.DEG <- decideTestsDGE(res_fls2vsKit_germfree, p.value=0.05, adjust.method = "BH", lfc=log2(2))
summary(is.DEG)
fls2_germfree_ID<-list(fls2_germfree_up     = rownames(is.DEG[is.DEG==1,1]),
                       fls2_germfree_down   = rownames(is.DEG[is.DEG==-1,1]))
p<-glimmaVolcano(res_fls2vsKit_germfree);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "fls2_vs_Kit_germfree.html"))


fls2vsKit_microbiota <- makeContrasts(fls2_microbiota-kit_microbiota, levels=design)
res_fls2vsKit_microbiota <- glmQLFTest(fit, contrast=fls2vsKit_microbiota )
is.DEG <- decideTestsDGE(res_fls2vsKit_microbiota, p.value=0.05, adjust.method = "BH", lfc=log2(2))
summary(is.DEG)
fls2_microbiota_ID<-list(fls2_microbiota_up   = rownames(is.DEG[is.DEG==1,1]),
                         fls2_microbiota_down = rownames(is.DEG[is.DEG==-1,1]))

fls2_DEG_ID<-c(fls2_germfree_ID,fls2_microbiota_ID);

####

############################# hupy edit (start) #######################
#  design <- model.matrix(~0+Group)
library(tidyr)
my.contrasts <- makeContrasts(
  Kit.MvsG = kit_microbiota - kit_germfree, # Genotype: WT Kitaake
  fls2.MvsG = fls2_microbiota - fls2_germfree, # Genotype: osfls2
  cerk.MvsG = cerk1_microbiota- cerk1_germfree, #Genotype: cerk1
  fls2vsKit_germfree =fls2_germfree -kit_germfree, # fls2 vs Kit, condition: germfree
  cerkvsKit_germfree = cerk1_germfree - kit_germfree, # cerk1 vs Kit, condition: germfree
  fls2vsKit_microbiota = fls2_microbiota - kit_microbiota, # fls2 vs Kit, condition: microbiota
  cerkvsKit_microbiota = cerk1_microbiota - kit_microbiota, # cerk1 vs Kit, condition: microbiota
  con_fls_x_microbiota = (fls2_microbiota - fls2_germfree) - (kit_microbiota - kit_germfree), # GLM model testing genes that differentially response to microbiota in fls2 and Kit (G X M)
  con_cerk_x_microbiota = (cerk1_microbiota- cerk1_germfree) - (kit_microbiota - kit_germfree), # GLM model testing genes that differentially response to microbiota in cerk1 and Kit (G X M)
  levels= design
) 

degnumlist<- data.frame(Var1=c(), Var2=c(), Freq=c(), group=c())
for(i in colnames(my.contrasts)){
  qlf <- glmQLFTest(fit, contrast=my.contrasts[,i])
  topTags(qlf, sort.by = "logFC", p.value = 0.05)
  is.DEG <- decideTestsDGE(qlf, p.value=0.05, adjust.method="BH", lfc=log2(2) )
  summary(is.DEG) 
  summ<- data.frame(summary(is.DEG) ); summ$group<- i
  degnumlist<- rbind(degnumlist, summ)
}
degnumlist<- degnumlist %>% pivot_wider(names_from = "Var1", values_from = "Freq")

write.csv(degnumlist,paste0(OUT_DIR,"hupy.edgeR_contrasts_degnumstats.csv"))

degnumlist %>% filter(grepl("Kit|fls2", group))

# Total degs in Fig5b.



#### extract constrast and decideTestsDGE of each contrast
colnames(my.contrasts)
paste("res", colnames(my.contrasts),sep="_")

res_Kit.MvsG <- glmQLFTest(fit, contrast=my.contrasts[,1])
is.DEG_Kit.MvsG <- decideTestsDGE(res_Kit.MvsG, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_Kit.MvsG);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "1_Kit.MvsG.html"))

res_fls2.MvsG <- glmQLFTest(fit, contrast=my.contrasts[,2])
is.DEG_fls2.MvsG<- decideTestsDGE(res_fls2.MvsG, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_fls2.MvsG);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "2_fls2.MvsG.html"))

res_cerk.MvsG <- glmQLFTest(fit, contrast=my.contrasts[,3])
is.DEG_cerk.MvsG<- decideTestsDGE(res_cerk.MvsG, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_cerk.MvsG);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "3_cerk.MvsG.html"))

res_fls2vsKit_germfree <- glmQLFTest(fit, contrast=my.contrasts[,4])
is.DEG_fls2vsKit_germfree<- decideTestsDGE(res_fls2vsKit_germfree, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_fls2vsKit_germfree);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "4_fls2vsKit_germfree.html"))

res_cerkvsKit_germfree <- glmQLFTest(fit, contrast=my.contrasts[,5])
is.DEG_cerkvsKit_germfree<- decideTestsDGE(res_cerkvsKit_germfree, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_cerkvsKit_germfree);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "5_cerkvsKit_germfree.html"))

res_fls2vsKit_microbiota <- glmQLFTest(fit, contrast=my.contrasts[,6])
is.DEG_fls2vsKit_microbiota<- decideTestsDGE(res_fls2vsKit_microbiota, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_fls2vsKit_microbiota);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "6_fls2vsKit_microbiota.html"))

res_cerkvsKit_microbiota <- glmQLFTest(fit, contrast=my.contrasts[,7])
is.DEG_cerkvsKit_microbiota<- decideTestsDGE(res_cerkvsKit_microbiota, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_cerkvsKit_microbiota);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "7_cerkvsKit_microbiota.html"))

res_con_fls_x_microbiota <- glmQLFTest(fit, contrast=my.contrasts[,8])
is.DEG_con_fls_x_microbiota<- decideTestsDGE(res_con_fls_x_microbiota, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_con_fls_x_microbiota);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "8_con_fls_x_microbiota.html"))

res_con_cerk_x_microbiota <- glmQLFTest(fit, contrast=my.contrasts[,9])
is.DEG_con_cerk_x_microbiota<- decideTestsDGE(res_con_cerk_x_microbiota, p.value=0.05, adjust.method="BH", lfc=log2(2) )
p<-glimmaVolcano(res_con_cerk_x_microbiota);
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "9_con_cerk_x_microbiota.html"))


#### prepare DEG list for Venn #######

########## osfls2_deg_list ########## 
is.DEG_fls2vsKit_microbiota
is.DEG_fls2vsKit_germfree

fls2_germfree_ID<-list(fls2_germfree_up   = rownames(is.DEG_fls2vsKit_germfree[is.DEG_fls2vsKit_germfree==1,1]),
                       fls2_germfree_down = rownames(is.DEG_fls2vsKit_germfree[is.DEG_fls2vsKit_germfree==-1,1]))

fls2_microbiota_ID<- list(fls2_microbiota_up   = rownames(is.DEG_fls2vsKit_microbiota[is.DEG_fls2vsKit_microbiota == 1 ,1]),
                          fls2_microbiota_down = rownames(is.DEG_fls2vsKit_microbiota[is.DEG_fls2vsKit_microbiota ==-1 ,1]))

fls2_x_microbiota_ID<-list(fls2_x_m_up   = rownames(is.DEG_con_fls_x_microbiota[is.DEG_con_fls_x_microbiota==1,1]),
                           fls2_x_m_down = rownames(is.DEG_con_fls_x_microbiota[is.DEG_con_fls_x_microbiota==-1,1]))

fls2_DEG_ID<-c(fls2_germfree_ID,fls2_microbiota_ID,fls2_x_microbiota_ID);
fls2_DEG_ID[c(1,3)]

fls2_DEG_allup<- c(fls2_germfree_ID$fls2_germfree_up, fls2_microbiota_ID$fls2_microbiota_up, fls2_x_microbiota_ID$fls2_x_m_up) #1317
fls2_DEG_alldown<- c(fls2_germfree_ID$fls2_germfree_down, fls2_microbiota_ID$fls2_microbiota_down, fls2_x_microbiota_ID$fls2_x_m_down)

num1<- fls2_DEG_allup %>% length()
num2<- fls2_DEG_allup %>% unique() %>% length()
num3<- fls2_DEG_alldown %>% length()
num4<- fls2_DEG_alldown %>% unique() %>% length()
num5<- c(fls2_DEG_allup %>% unique(), fls2_DEG_alldown %>% unique() )%>% length()
num6<- c(fls2_DEG_allup %>% unique(), fls2_DEG_alldown %>% unique() )%>% unique() %>% length()

mutant_all_degs<- data.frame(mutant=c(), up_dup=c(), up_derep=c(), down_dup=c(), down_derep=c(), all_dup=c(), all_derep=c())
mutant_all_degs<- rbind(mutant_all_degs, data.frame(mutant="fls2", up_dup=num1, up_derep=num2, down_dup=num3, down_derep=num4, all_dup=num5, all_derep=num6))

########## cerk1_deg_list ##########

cerk1_germfree_ID<-list(cerk1_germfree_up   = rownames(is.DEG_cerkvsKit_germfree[is.DEG_cerkvsKit_germfree==1,1]),
                       cerk1_germfree_down = rownames(is.DEG_cerkvsKit_germfree[is.DEG_cerkvsKit_germfree==-1,1]))

cerk1_microbiota_ID<- list(cerk1_microbiota_up   = rownames(is.DEG_cerkvsKit_microbiota[is.DEG_cerkvsKit_microbiota == 1 ,1]),
                          cerk1_microbiota_down = rownames(is.DEG_cerkvsKit_microbiota[is.DEG_cerkvsKit_microbiota ==-1 ,1]))

cerk1_x_microbiota_ID<-list(cerk1_x_m_up   = rownames(is.DEG_con_cerk_x_microbiota[is.DEG_con_cerk_x_microbiota==1,1]),
                            cerk1_x_m_down = rownames(is.DEG_con_cerk_x_microbiota[is.DEG_con_cerk_x_microbiota==-1,1]))

cerk1_DEG_ID<-c(cerk1_germfree_ID,cerk1_microbiota_ID,cerk1_x_microbiota_ID);
cerk1_DEG_ID[c(1,3)]

cerk1_DEG_allup<- c(cerk1_germfree_ID$cerk1_germfree_up, cerk1_microbiota_ID$cerk1_microbiota_up, cerk1_x_microbiota_ID$cerk1_x_m_up) #1317
cerk1_DEG_alldown<- c(cerk1_germfree_ID$cerk1_germfree_down, cerk1_microbiota_ID$cerk1_microbiota_down, cerk1_x_microbiota_ID$cerk1_x_m_dowm)

num1<- cerk1_DEG_allup %>% length()
num2<- cerk1_DEG_allup %>% unique() %>% length()
num3<- cerk1_DEG_alldown %>% length()
num4<- cerk1_DEG_alldown %>% unique() %>% length()
num5<- c(cerk1_DEG_allup %>% unique(), cerk1_DEG_alldown %>% unique() )%>% length()
num6<- c(cerk1_DEG_allup %>% unique(), cerk1_DEG_alldown %>% unique() )%>% unique() %>% length()

mutant_all_degs<- rbind(mutant_all_degs, data.frame(mutant="cerk1", up_dup=num1, up_derep=num2, down_dup=num3, down_derep=num4, all_dup=num5, all_derep=num6))
degnumlist
write.csv(mutant_all_degs,paste(OUT_DIR,"hupy.mutant_all_degs.csv",sep=""), quote=F)

######################## hupy edit (end) #########################

#### Draw Venn diagram for fls2-associated DEGs ######

# down-regulated genes
p1<-ggvenn(
  fls2_DEG_ID[c(2,4,5,6)],
  fill_color = c("#0073C2FF", "#CD534CFF", "#CD534CFF", "#868686FF"),
  show_percentage = TRUE, auto_scale = FALSE,
  stroke_size = 0.5, set_name_size = 4
)

#up-regulated genes
p2<-ggvenn(
  fls2_DEG_ID[c(1,3,5,6)],
  fill_color = c("#0073C2FF", "#CD534CFF", "#CD534CFF", "#868686FF"),
  show_percentage = TRUE, auto_scale = FALSE,
  stroke_size = 0.5, set_name_size = 4
)

p<-plot_grid(p1,p2)
p
ggsave(file.path(OUT_DIR,"venn_osfls2.pdf"),p)

UpSetR::fromList(fls2_DEG_ID)
pdf(file.path(OUT_DIR, "fls2_UpSet_draw.pdf"))
upset(fromList((fls2_DEG_ID)), keep.order =T)
dev.off()

### Heatmaps #########
length(unique(sort(unlist(fls2_DEG_ID))))  # 2090 genes associated to fls2
fls2_DEG_cpm<-cpms[unique(sort(unlist(fls2_DEG_ID))),]
p<-pheatmap(fls2_DEG_cpm[,1:12],show_rownames = F, scale="row",main = "CPM of DEGs associated to FLS2");

pdf(file=file.path(OUT_DIR,"fls2_cpm_heatmap.pdf"),width=4,height = 10)
p
dev.off()
## Combined the logFC and pvalue together
## Check the row names in different datasets
table(rownames(res_fls2vsKit_germfree$table)==rownames(res_fls2vsKit_microbiota$table));

#Extract logFC values to a dataframe and a file
GOI<- read.csv("out_edgeR2/Figure5ht_fls2_stats_GOI.csv", header = T, row.names = 1) %>%
  select(starts_with("GOI"))
unique(GOI$GOI_PR)
unique(GOI$GOI_Auxin)
unique(GOI$GOI_WRKY)
table(rownames(GOI) == rownames(fls2_stat))

logfc_1<- res_fls2vsKit_germfree$table
logfc_2<- res_fls2vsKit_microbiota$table
logfc_3<- res_fls_x_microbiota$table
logfc_4<- res_Kit.MvsG$table
cpm_zscore<- cpms %>% t() %>% scale() %>% t() %>%as.data.frame() %>%
  select(starts_with("fls2"), starts_with("kit"))
names(cpm_zscore)<- paste(names(cpm_zscore), "_cpm_zscore", sep="")

colnames(logfc_1)<- paste0("fls2_G_", colnames(logfc_1))
colnames(logfc_2)<- paste0("fls2_M_", colnames(logfc_2))
colnames(logfc_3)<- paste0("fls2_GxM_", colnames(logfc_3))
colnames(logfc_4)<- paste0("Kit_", colnames(logfc_4))

fls2_stat<-cbind(logfc_3,logfc_2,logfc_1,logfc_4,cpm_zscore)[unique(sort(unlist(fls2_DEG_ID))),] %>%
  arrange(fls2_GxM_logFC)#2090 degs

fls2_stat %>% 
  select(ends_with("logFC")) %>%summary
fls2_stat %>% 
  select(ends_with("_cpm_zscore")) %>%View()
col_logFC<- colorRamp2(c(-12,-6,-2,-1,0,1,2,6,12), colorRampPalette(rev(brewer.pal(n=10, name="BrBG")))(9))
col_logCPM<- colorRamp2(c(-12,-6,-2,-1,0,1,2,6,12), colorRampPalette(rev(brewer.pal(n=10, name="RdBu")))(9))
col_GOI = colorRamp2(c(1,0), c("blue","white"))

ht1<- fls2_stat %>% 
  select(ends_with("logFC")) %>% 
  as.matrix() %>%
  Heatmap(name= "logFC",cluster_rows=F, cluster_columns = F,
          col= col_logFC)

ht2<- fls2_stat %>% 
  select(ends_with("_cpm_zscore")) %>%
  as.matrix() %>%
  Heatmap(name= "CPM_zscore",cluster_rows=F, cluster_columns = F,
          col= col_logCPM)

ht3<- as.matrix(GOI) %>% 
  Heatmap(name= "GOI",cluster_rows=F, cluster_columns = F,col=col_GOI)

pdf(paste0(OUT_DIR,"Fig5.ht_fls2_degs.pdf"))
ht1+ht2+ht3
dev.off()        

######################## 20240510 edit END ######################
color<-c(rev(brewer.pal(4, "Blues")),"#FFFFFF","#FFFFFF",brewer.pal(4,"Reds"))
#color<-c(brewer.pal(20, "Blues"),"#000000","#000000",rev(brewer.pal(4,"Reds")))  # REVERSE THE COLOR
pheatmap(fls2_logFC,show_rownames = F, scale="none",cluster_cols=T,color = color,
         breaks=c(-8,-6,-4,-2,-1,0,1,2,4,6,8));

# order the DEGs by logFC
fls2_logFC<-fls2_logFC %>% arrange(desc(GXM), Germfree,Microbiota)
#annotation_row = Os_entrez[rownames(fls2_logFC),c(1,18)]
p<-pheatmap(fls2_logFC,show_rownames = F, 
         scale="none",cluster_cols=T,cluster_rows=F,color = color,
         breaks=c(-8,-6,-4,-2,-1,0,1,2,4,6,8),
         legend_breaks = c(-8,-6,-4,-2,-1,0,1,2,4,6,8),
         #annotation_row = annotation_row,
         main="DEGs in osfls2 mutant")
pdf(file=file.path(OUT_DIR,"fls2_logFC.pdf"),width=4,height = 10)
p
dev.off()

## absent/present heatmap(+1,-1, 0)
d<-cbind(decideTestsDGE(res_fls2vsKit_microbiota, p.value=0.05, adjust.method="BH", lfc=log2(2)),
      decideTestsDGE(res_fls2vsKit_germfree, p.value=0.05, adjust.method="BH", lfc=log2(2) ),
      decideTestsDGE(res_fls_x_microbiota, p.value=0.05, adjust.method="BH", lfc=log2(2)))

head(d)
colnames(d)<-c("microbiota","germfree","GxM")
d<-d[unique(sort(unlist(fls2_DEG_ID))),]
pheatmap(d,show_rownames = F, scale="none",cluster_cols=T,color = c("red","grey","green"));

#### 20240530 BiNGO_GO point_fig ######
library(dplyr)
library(stringr)

fls2_bingo_res<- readxl::read_excel("./BiNGO/fls2_vs_Kit_integrated_GO.xlsx")
fls2_bingo_res$GO <- paste("GO ",
                           str_pad(fls2_bingo_res$GO_ID, width = 8, side = "left", pad = "0"),
                           ": ",
                           fls2_bingo_res$`GO Description`)
pattern <- "protein|radiation|nucleobase|transcription|RNA|macromolecule|histone|mono|expression|biosynthetic"
plotdata <- fls2_bingo_res %>% 
  #filter(GeneNo >= 0) %>%
  filter(!str_detect(GO, pattern))
count_data <- plotdata %>% count(GO, name = "n")
merged_data <- plotdata %>% inner_join(count_data, by = "GO")

p1<- merged_data %>% 
  select(Group, GeneNo, `GO Description`, Corrected_pval, n) %>% 
  filter(GeneNo > 15 ) %>%
  filter(Group %in% c("fls2_germfree_1303", "fls2_microbiota_514", "fls2_x_microbiota_1181")) %>%
  mutate(Group=factor(Group, levels=c("fls2_x_microbiota_1181","fls2_microbiota_514","fls2_germfree_1303"))) %>%
  ggplot(aes(x=Group, y=`GO Description`))+
  geom_point(aes(size=GeneNo, color=Corrected_pval))+
  theme_bw()+
  scale_color_gradient(high = "#D8BFD8", low = "#6C2E53")+
  theme(axis.text = element_text(color="black",size=12),
        axis.text.x = element_text(angle = 90, hjust=1))+
  labs(title = "GeneNo>15",x="", y="GO term")

p2<- merged_data %>% 
  select(Group, GeneNo, `GO Description`, Corrected_pval, n) %>% 
  filter(GeneNo > 15 ) %>%
  filter(!(Group %in% c("fls2_germfree_1303", "fls2_microbiota_514", "fls2_x_microbiota_1181"))) %>%
  ggplot(aes(x=Group, y=`GO Description`))+
  geom_point(aes(size=GeneNo, color=Corrected_pval))+
  theme_bw()+
  scale_color_gradient(high = "#D8BFD8", low = "#6C2E53")+
  theme(axis.text = element_text(color="black",size=12),
        axis.text.x = element_text(angle = 30, hjust=1))+
  scale_x_discrete(expand = c(0.3, 0.3)) + 
  labs(title = "GeneNo>15",x="", y="GO term")

library(cowplot)
combined_plot <- plot_grid(p1 , p2 ,ncol = 1, align = "v", rel_heights  = c(5, 6))
print(combined_plot)
ggsave("./Fig5.BiNGO_point_GOn_thred3.pdf", width=7, height=9)


ggsave("tmpfig.pdf", width=18,height=8)



########

#########################################################################
## Output the results
cbind(res_fls2vsKit_germfree$table, res_fls2vsKit_microbiota$table,res_fls_x_microbiota$table)[unique(sort(unlist(fls2_DEG_ID))),]
fls2_all_contrast<-makeContrasts(fls2DEG_germfree   = fls2_germfree-kit_germfree,
              fls2DEG_microbiota = fls2_microbiota-kit_microbiota,
              fls2DEG_x_microbiota = (fls2_microbiota-fls2_germfree)-(kit_microbiota-kit_germfree),
              levels=design)
res_fls2_all<-glmQLFTest(fit, contrast = fls2_all_contrast);
topTags(res_fls2_all);
is.DEG<-decideTestsDGE(res_fls2_all,p.value=0.05, adjust.method = "BH",lfc = log2(2))
summary(is.DEG)     #2015 genes
## get fls2 DEGs
length(unique(sort(unlist(fls2_DEG_ID))))  # 2090 genes

fls2_out_df <-data.frame(topTags(res_fls2_all, n=30000)[unique(sort(unlist(fls2_DEG_ID))),]) 
write.table(fls2_out_df, file = file.path(OUT_DIR, "fls2.DEG.txt"))
############################################################################################
##
##  processing cerk1 RNA-seq
##

cerk1vsKit_germfree <- makeContrasts(groupcerk1_germfree-groupkit_germfree, levels=design)
res_cerk1vsKit_germfree <-glmQLFTest(fit,contrast = cerk1vsKit_germfree)
topTags(res_cerk1vsKit_germfree,sort.by="logFC")
is.DEG<-decideTestsDGE(res_cerk1vsKit_germfree,p.value=0.05, adjust.method = "BH",lfc = log2(2))
summary(is.DEG)

cerk1_germfree_ID<-list(cerk1_germfree_up = row.names(is.DEG[is.DEG==1,1]),
                        cerk1_germfree_down= row.names(is.DEG[is.DEG==-1,1]))

p<-glimmaVolcano(res_cerk1vsKit_germfree);
p
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "cerk1_vs_Kit_germfree.html"))

## cerk1 DEGs in microbiota condition
cerk1vsKit_microbiota <- makeContrasts(groupcerk1_microbiota-groupkit_microbiota, levels=design)
res_cerk1vsKit_microbiota <- glmQLFTest(fit, contrast = cerk1vsKit_microbiota)
topTags(res_cerk1vsKit_microbiota);
is.DEG<-decideTestsDGE(res_cerk1vsKit_microbiota,p.value=0.05, adjust.method = "BH",lfc=log2(2))
summary(is.DEG)
cerk1_microbiota_ID<- list(cerk1_microbiota_up =row.names(is.DEG[is.DEG==1,1]),
                           cerk1_microbiota_down = row.names(is.DEG[is.DEG==-1,1]))
p<-glimmaVolcano(res_cerk1vsKit_microbiota)
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "cerk1_vs_Kit_microbiota.html"))

# cerk1 DEGs in GXM conditions
con_cerk1_x_microbiota <- makeContrasts(
                           (groupcerk1_microbiota-groupcerk1_germfree)-(groupkit_microbiota-groupkit_germfree),
                            levels=design)
                           
res_cerk1_x_microbiota<-glmQLFTest(fit,contrast = con_cerk1_x_microbiota)
is.DEG<-decideTestsDGE(res_cerk1_x_microbiota,p.value=0.05, adjust.method = "BH", lfc = log2(2))
summary(is.DEG);
cerk1_x_M_ID<- list(cerk1_x_M_up =row.names(is.DEG[is.DEG==1,1]),
                    cerk1_x_M_down = row.names(is.DEG[is.DEG==-1,1]))

p<-glimmaVolcano(res_cerk1_x_microbiota)
htmlwidgets::saveWidget(p, file.path(OUT_DIR, "cerk1_x_microbiota.html"))

## Venn diagram for cerk1
cerk1_DEG_ID<-c(cerk1_germfree_ID, cerk1_microbiota_ID, cerk1_x_M_ID)
ggvenn(cerk1_DEG_ID[c(1,3,5,6)])
ggvenn(cerk1_DEG_ID[c(2,1,3,4)])
ggvenn(cerk1_DEG_ID[c(2,1,3,4)])
sapply(cerk1_DEG_ID, length); # counting the DEG number

##### Draw Venn diagram for cerk1-associated DEGs

# down-regulated genes
p1<-ggvenn(
  cerk1_DEG_ID[c(2,1,3,4)],
  fill_color = c("#0073C2FF", "#CD534CFF", "#CD534CFF", "#868686FF"),
  show_percentage = TRUE, auto_scale = FALSE,
  stroke_size = 0.5, set_name_size = 4
)
ggsave(file.path(OUT_DIR,"venn_cerk1.pdf"),p1)

UpSetR::fromList(cerk1_DEG_ID)
pdf(file.path(OUT_DIR, "cerk1_UpSet_draw.pdf"))
upset(fromList((cerk1_DEG_ID)), keep.order =T)
dev.off();

## Take DEGs in both fls2 and cerk1 as core-responsive genes
all_DEG<-list(FLS2_UNI_DEG=sort(unlist(fls2_DEG_ID)),
              CERK1_UNI_DEG=sort(unlist(cerk1_DEG_ID)))
ggvenn(all_DEG)
pdf(file.path(OUT_DIR, "cerk1_fls2_UpSet_draw.pdf"))
upset(fromList(c(fls2_DEG_ID[c(1,3,5)], cerk1_DEG_ID[c(1,3,5)])))
upset(fromList(c(fls2_DEG_ID[c(2,4,6)], cerk1_DEG_ID[c(2,4,6)])))
dev.off()




##########################################################################
## Function enrichment analysis
##  Not finish

## Access the NCBI entrez database
hub<-AnnotationHub()
query(hub,c("osa"))
Os_hub<-query(hub,c("Oryza sativa"));

## SYMBOL column contain LOCXXXX annotations
test.go<-enrichGO(all_DEG[[1]],
         OrgDb = Os_hub[["AH114588"]],
         keyType = "SYMBOL",
         ont = "ALL")
head(test.go)
dim(test.go)
dotplot(test.go)

test.go<-enrichGO(fls2_DEG_ID[[3]],
                  OrgDb = Os_hub[["AH114588"]],
                  keyType = "SYMBOL",
                  universe = substr(tx2gene$V2,4,20),
                  ont = "ALL")
dotplot(test.go)
data.frame(test.go)

test.go<-enrichGO(unlist(fls2_germfree_ID),
                  OrgDb = Os_hub[["AH114588"]],
                  keyType = "SYMBOL",
                  ont = "ALL")
dotplot(test.go)
cnetplot(test.go, fls2_logFC$Germfree)


test.go<-gseGO(sort(as.integer(test),decreasing = T),
      OrgDb = Os_hub[["AH114588"]],
      maxGSSize = 3000)
head(test.go)
goplot(test.go,howCategory = 5)

#################################
## Using downloaded GO annotation for enrichment analysis
Os_go<-read.table(file="~/project/pubdata/NIP_ncbi_dataset/data/GCF_001433935.1/GCF_001433935.1_IRGSP-1.0_gene_ontology.gaf", head=F, skip = 9)
Os_go<-Os_go[,c(3,5)]

############################33
##
search_kegg_organism('osa', by='kegg_code')

test_kegg <- enrichKEGG(gene  = substr(all_DEG[[1]],4,20),
                 organism     = 'osa',
                 universe = substr(tx2gene$V2,4,20),
                 maxGSSize = 3000,
                 pvalueCutoff = 0.05)

data.frame(test_kegg)

length(test)
test.ID<-paste("LOC",unlist(strsplit(data.frame(test_kegg)[3,"geneID"], "/")), sep="")
pheatmap(fls2_DEG_cpm[test.ID,1:12],show_rownames = T, scale="none",cluster_cols=F,cluster_rows=F,color = color);

color<-c(rev(brewer.pal(4, "Blues")),"#FFFFFF","#FFFFFF",brewer.pal(4,"Reds"))
#color<-c(brewer.pal(20, "Blues"),"#000000","#000000",rev(brewer.pal(4,"Reds")))  # REVERSE THE COLOR
## genes associated to a enriched KEGG pathway
annotation_row = Os_entrez[test.ID,c(18)]
kegg_enriched<-fls2_logFC[test.ID,];
rownames(kegg_enriched)[annotation_row!="-"]<-annotation_row[annotation_row!="-"]
pheatmap(kegg_enriched,show_rownames = T, scale="none",cluster_cols=F,color = color,
         breaks=c(-8,-6,-4,-2,-1,0,1,2,4,6,8));

data.frame(test_kegg)
Os_entrez[test.ID,]

library(ggkegg)
graph <- ggkegg::pathway("osa00904", use_cache=TRUE)
graph


