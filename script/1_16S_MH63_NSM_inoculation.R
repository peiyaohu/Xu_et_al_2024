#### Script1 #####
# Obj: 16S of MH63_NSM_inoculation, including:
# > ExtFig2A-C. diversity
# > Fig2A. heatmap


### 0. Set path and load data ######
library(reshape2)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(vegan)
library(tidyverse)
library(agricolae)

path<- "~/Documents/Research/LabMember/xjw/Publish/DataAvailibility/Code/"
setwd(path)
output.path<- "/Users/hupeiyao/Documents/Research/LabMember/xjw/Publish/DataAvailibility/Code/"
source("save_pheatmap.R")

main_theme =   theme(axis.line.x=element_line(size=.5, colour="black"),
                     axis.line.y=element_line(size=.5, colour="black"),
                     axis.ticks=element_line(color="black"),
                     axis.text=element_text(color="black", size=12))


asv_tax<- read.csv("source/1_16S_MH63_NSM/asv_tax.csv",header = T,row.names = 1)
asv_tab<- read.csv("source/1_16S_MH63_NSM/asv_tab.csv",header = T,row.names = 1)
metadata<- read.table("source/1_16S_MH63_NSM/sample-metadata.tsv",header = T,row.names = 1,sep="\t")

### 1. construct phyloseq object; filt low-quality reads and rarefy #######
ps<- phyloseq(otu_table(asv_tab,taxa_are_rows = T),
              tax_table(as.matrix(asv_tax)),
              sample_data(metadata))  #7312
### 1) filt low-quality
# filt1: remove host contamination
idx1<- grep("Mito",asv_tax$Family);idx1 #23
idx2<- which(asv_tax$Kingdom!="Bacteria");idx2 #6
idx3<- which(is.na(asv_tax$Phylum));idx3 #130

idx<- c(grep("Mito",asv_tax$Family),
        which(asv_tax$Kingdom!="Bacteria"), #6
        which(is.na(asv_tax$Phylum))) #130
#sum up: 159
length(idx)

ps.1<- prune_taxa(!rownames(asv_tax) %in% rownames(asv_tax)[idx1],ps);ps.1 #7289/7312
ps.2<- prune_taxa(!rownames(otu_table(ps.1)) %in% rownames(asv_tax)[idx2],ps.1);ps.2 #7283/7289
ps.3<- prune_taxa(!rownames(otu_table(ps.2)) %in% rownames(asv_tax)[idx3],ps.2);ps.3 #7153/7283

## statistic: reads of each step
stat<- data.frame(raw=colSums(asv_tab))
stat$removeMito<- sample_sums(ps.1)[rownames(stat)]
stat$removeNonBact<- sample_sums(ps.2)[rownames(stat)]
stat$removePhylumNA<- sample_sums(ps.3)[rownames(stat)]
head(stat)

# filt1: remove all non-bacterial seqs(phylumNA; Mito; NonBact)
ps.filt<- prune_taxa(!rownames(asv_tax) %in% rownames(asv_tax)[idx],ps);ps.filt #7153
rm(ps.1,ps.2,ps.3,idx1,idx2,idx3,idx)
# filt2: rowSums<3
ps.filt<- prune_taxa(rowSums(otu_table(ps.filt))>3,ps.filt);ps.filt #6205/7153

### 2) rarefaction
set.seed(26)
# ascending order and set rarefaction depth
sample_sums(ps.filt)[order(sample_sums(ps.filt))]
#XF1_3 o1_xgl_soil_1 s1_xgl_soil_3         XF1_2 o1_xgl_soil_2 o1_xgl_soil_3 s1_xgl_soil_1 s1_xgl_soil_2 
#10224         17820         25900         31796         36333         40599         42998         43930 
#o2_xgl_soil_3 s2_xgl_soil_2 s1_xgl_root_1 s1_xgl_root_3 s2_xgl_soil_1 s1_xgl_root_2 o2_xgl_root_1 o1_xgl_root_2 
#45552         48289         54650         57960         60980         64097         64664         66418 
#XF2_1 o2_xgl_root_3 o1_xgl_root_3         XF1_1 s2_xgl_root_1         XF2_2 o1_xgl_root_1 o2_xgl_soil_2 
#66730         66772         68398         71730         75694         77034         86885         91002 
#s2_xgl_root_2         XF2_3 o2_xgl_soil_1 o2_xgl_root_2 s2_xgl_root_3 s2_xgl_soil_3 
#92799         98163        103169        105171        114802        144254 

### Final to output
ps.rarefy<- rarefy_even_depth(ps.filt,sample.size =31796 ) #5909/6205
# 3 samples were removed caused fewer reads than `sample.size`(=31796)
# o1_xgl_soil_1 s1_xgl_soil_3 XF1_3

asvtab<- data.frame(t(otu_table(ps.rarefy)));dim(asvtab)
asvtax<- data.frame(t(tax_table(ps.rarefy)));dim(asvtax)
asvtab_tax<- data.frame(t(asvtab));asvtab_tax<- cbind(asvtab_tax, asvtax)
write.csv(asvtab_tax,paste0(output.path,"Fig2.Table0.ps.rarefy.asvtab_5909.csv"))

#### ExtFig2A-C. Diversity #################
##### ExtFig2A-B ####
alpha<- estimate_richness(ps.rarefy)
alpha<- cbind(metadata,alpha[rownames(metadata),])
alpha$SoilType<- factor(alpha$SoilType,levels = c("Organic","Loam","Input"))
alpha$seqCompartment<- factor(alpha$seqCompartment,levels = c("Input","Rhizo","Root"))
alpha$groups<- factor(alpha$groups,levels = c("Input" ,"LoamRhizo" ,"LoamRoot",
                                              "OrganicRhizo", "OrganicRoot"))

# 1- Group by Compartment
## Combine 2 bacthes together.
m<- c("Shannon","Observed","Chao1")
for(n in 1:3){
  i=m[n]
  model<- aov(alpha[[i]]~groups,data=alpha)
  Tukey_HSD<-TukeyHSD(model,ordered = T,conf.level = .95)
  names(Tukey_HSD)
  Tukey_HSD_table<-as.data.frame(Tukey_HSD$groups)
  
  write.table(Tukey_HSD_table,file=paste0(output.path,"ExtFig2A.",i,"-test_2batches.tsv"),sep="\t")

  library(agricolae)
  library(magrittr)
  out<-LSD.test(model,"groups",p.adj="none")
  head(out);names(out)
  stat<-out$groups
  
  alpha %<>% filter(!is.na(alpha$Observed)) 
  alpha$stat<-stat[as.character(alpha$groups),]$groups
  max<-max(alpha[,i])
  min<-min(alpha[,i])
  library(dplyr)
  x<-alpha[,c("groups",i)]
  y <- x %>% 
    group_by(groups) %>%
    summarize_(Max=paste("max(",i,")", sep=""))
  y<-as.data.frame(y)
  rownames(y)<-y$groups
  alpha$y<-y[as.character(alpha$groups),]$Max+(max-min)*0.05
  
  p<-ggplot(alpha,aes(x=groups,y=get(i),color=groups))+
    geom_boxplot(alpha=1,outlier.size = 1,size=.7,alpha=.7,width=0.5,fill="transparent")+
    labs(x="Groups",y=paste(i,"index"))+theme_classic()+
    geom_text(data=alpha,aes(x=groups,y=y,color=groups,label=stat),size=5)+
    geom_jitter(position=position_jitter(.17),size=1,alpha=.7)+
    theme(axis.line.x=element_line(size=.5, colour="black"),
          axis.line.y=element_line(size=.5, colour="black"),
          axis.ticks=element_line(color="black"),
          axis.text=element_text(color="black", size=12))
  
  color=c("#1F4E79", "#9E7C61","#00BFC4","#9E7C61","#00BFC4")
  p<- p+scale_color_manual(values=color)
  ggsave(paste0(output.path,"ExtFig2A.alpha_",i,"_2batches.pdf",sep=""),p,width=5,height=3)
  write.table(alpha,paste0(output.path,paste0("ExtFig2A.",i,"-plot_2bacthes.tsv")),sep="\t")
}

##### ExtFig2C Beta diversity #########
bray_all<- vegdist(asvtab,method = "bray")

# PCoA to plot:
bray=bray_all
submeta<- metadata[colnames(as.matrix(bray_all)),]
perm1<-  adonis2(bray_all ~ seqCompartment + SoilType + Batch,data=submeta, permutations = 9999)
perm2<-  adonis2(bray_all ~   SoilType +seqCompartment+ Batch,data=submeta, permutations = 9999)

df<- data.frame(perm1)
df<- rbind(df,data.frame(perm2))
write.csv(df,paste0(output.path,"ExtFig2C.Tab_PERMANOVA.csv",sep=""))

#### NMDS
df_nmds <- metaMDS(bray_all, k = 2)
# stress、points & species
summary(df_nmds)
df_nmds_stress <- df_nmds$stress
df_nmds_stress
stressplot(df_nmds)

df_points <- as.data.frame(df_nmds$points)
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')

# add group info
df_points$seqCompartment<- factor(metadata[rownames(df_points),"seqCompartment"],levels=c("Input","Rhizo","Root"))
df_points$SoilType=metadata[rownames(df_points),"SoilType"]

p=ggplot(df_points,aes(x=NMDS1,y=NMDS2,color=seqCompartment,shape=SoilType))+
  geom_point(alpha=.7,size=3)+
  labs(x=paste("NMDS1 (",format(100 *eig[1]/sum(eig),digits = 4),"%)",sep=""),
       y=paste("NMDS2 (",format(100 *eig[2]/sum(eig),digits = 4),"%)",sep=""),
       title=paste("bray_curtis"," NMDS",sep=""))+
  scale_color_manual(values=color)+
  theme_classic()+main_theme+
  xlim(c(-2,2.5))+ylim(c(-2,2))
ggsave(paste0(output.path,"ExtFig2C.NMDS_all_legend.pdf",sep=""),width = 5,height = 3)
write.csv(df_points,paste0(output.path,"ExtFig2C.Bray_NMDS.csv",sep=""))


#### Figure 2A. Heatmap #######
phylum_color<- c("#ac8daf","#ddb6c6","#f1d4d4","#85A3BF","#f19584","#4F9153","#6C5B7B")

###### 1) WGCNA: modularity ==> preparation for heatmap's row_cluster
library(WGCNA)
library(reshape2)
library(viridis)
library(RColorBrewer)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(asvtab, powerVector = powers, verbose = 5)
# set 8
net = blockwiseModules(asvtab, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "xjwCup16TOM",
                       verbose = 3)
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FxjwCup16S-02-networkConstruction-auto.RData")

## Modules:
module <- data.frame(net$colors) 
colnames(module)<- "Modules"
library(dplyr)
module<- cbind(module,asv_tax[rownames(module),])


ps.rarefy.sandy.root<- prune_samples(sample_data(ps.rarefy)$SoilType=="Loam" &sample_data(ps.rarefy)$seqCompartment=="Root" ,ps.rarefy)
ps.rarefy.sandy.root<- prune_taxa(rowSums(otu_table(ps.rarefy.sandy.root))>0,ps.rarefy.sandy.root) #604
module$rarefy_sandyroot<- ifelse(rownames(module) %in% taxa_names(ps.rarefy.sandy.root),1,0)

ps.rarefy.sandy.soil<- prune_samples(sample_data(ps.rarefy)$SoilType=="Loam" &sample_data(ps.rarefy)$seqCompartment=="Rhizo" ,ps.rarefy)
ps.rarefy.sandy.soil<- prune_taxa(rowSums(otu_table(ps.rarefy.sandy.soil))>0,ps.rarefy.sandy.soil) #1539
module$rarefy_sandysoil<- ifelse(rownames(module) %in% taxa_names(ps.rarefy.sandy.soil),1,0)

ps.rarefy.sandy.input<- prune_samples(sample_data(ps.rarefy)$seqCompartment=="Input",ps.rarefy)
ps.rarefy.sandy.input<- prune_taxa(rowSums(otu_table(ps.rarefy.sandy.input))>0,ps.rarefy.sandy.input) #2700
module$rarefy_input<- ifelse(rownames(module) %in% taxa_names(ps.rarefy.sandy.input),1,0)


module_stat<- module[,c(1,8:10)] %>% group_by(Modules) %>% summarise_each(sum)%>% data.frame()
#write.csv(module_stat,"4.WGCNA/Tab1_module_stats_sandy.csv")


stat<- module_stat %>% select(rarefy_sandyroot,rarefy_sandysoil,rarefy_input) %>% apply(2,function(x)x/sum(x))
rownames(stat)<- paste0("Module",module_stat[,1])
module_stat_melt<- melt(stat)
Levels=c("rarefy_input","rarefy_sandysoil","rarefy_sandyroot")
module_stat_melt$Var2<- factor(module_stat_melt$Var2,levels =Levels)


######## 2) Arrange ASVs by Modules
module_arrange<- module %>% arrange(Modules) 
asv_module<- data.frame(t(asvtab))[rownames(module_arrange),]
samples<- rownames(asvtab)
sample_order<- c(samples[grep("XF",samples)],
                 samples[grep("soil",samples)],
                 samples[grep("root",samples)])
dim(asv_module)

asv_module<- data.frame(t(asvtab))[rownames(module_arrange),]
samples<- rownames(asvtab)
sample_order<- c(samples[grep("XF",samples)],
                 samples[grep("soil",samples)],
                 samples[grep("root",samples)])
dim(asv_module)


# sandy, organic individually.
asv_module<- data.frame(t(asvtab))
samples<- rownames(asvtab)
sample_order<- c(samples[grep("XF",samples)],
                 samples[grep("soil",samples)],
                 samples[grep("root",samples)])

# heatmap's value
asv_module_mean<- data.frame(ASV=rownames(asv_module))
compartment<- c("XF","soil","root")
for(i in compartment){
  if(i == "soil"){
    meandf<- asv_module %>% select(matches(i)) 
    meandf_o<- meandf[,1:5] %>% apply(1,mean) %>% as.data.frame()
    meandf_s<- meandf[,6:10] %>% apply(1,mean) %>% as.data.frame()
    meandf<- cbind(meandf_o,meandf_s)
    colnames(meandf)<- c(paste("organic",i,"mean",sep="_"),paste("loam",i,"mean",sep="_"))
  }else if(i =="root"){
    meandf<- asv_module %>% select(matches(i)) 
    meandf_o<- meandf[,1:6] %>% apply(1,mean) %>% as.data.frame()
    meandf_s<- meandf[,7:12] %>% apply(1,mean) %>% as.data.frame()
    meandf<- cbind(meandf_o,meandf_s)
    colnames(meandf)<- c(paste("organic",i,"mean",sep="_"),paste("loam",i,"mean",sep="_"))
    
    }else{
    meandf<- asv_module %>% select(matches(i)) %>% apply(1,mean) %>% as.data.frame()
    colnames(meandf)<- paste(i,"mean",sep="_")
  }
  asv_module_mean<- cbind(asv_module_mean,meandf)
}
asv_module_mean<- asv_module_mean[,-1]

# Add taxonomy
tax_data<- asv_tax[rownames(asv_module_mean),c("Phylum","Class","Family")]
tax_data$Phylum[is.na(tax_data$Phylum)]="Unknown"
idx<- tax_data$Phylum=="Proteobacteria"
## choose top7 phyla to display
tax_data[idx,"Phylum"]=tax_data[idx,"Class"]
asv_module_mean$Phylum<- tax_data[rownames(asv_module_mean),"Phylum"]

## Generating new variable `tmp`; Add `Phylum` column.
tmp<- asv_module_mean %>% group_by(Phylum) %>% summarise_each(sum) %>% data.frame()
tmp$Phylum[is.na(tmp$Phylum)]="Unknown"
rownames(tmp)<- tmp$Phylum;tmp<- data.frame(tmp[,-1])
tmp$rowsum<- rowSums(tmp);tmp<- tmp %>% arrange(rowsum)

## Calculation based on tmp: each phyla RA_percentage(%):
## : rowSum(XF + o_s + s_s + o_r + s_r)
sum(tmp$rowsum)/5
phyla_RA_percentage<- data.frame(phyla=rownames(tmp),
                                 percentage=tmp$rowsum/(5*10276)*100);
phyla_RA_percentage<- phyla_RA_percentage[ order(phyla_RA_percentage$percentage,decreasing = T),]
rownames(phyla_RA_percentage)<- phyla_RA_percentage$phyla
top9<- rownames(phyla_RA_percentage)[1:9]
write.csv(phyla_RA_percentage,paste0(output.path,"Fig2A.Sup_phyla_RA_percentage.csv"))

# final choice
top9
# Group other phyla into `Others`
phyla<- rownames(tmp)
idx<- asv_module_mean$Phylum %in% top9
asv_module_mean$Phylum[!idx]<- "Others"
# Substitute all phyla
# Then for-loop:
top10<- c(top9,"Others")
# order top10:
proteo<- sort(top10[grep("proteo",top10)])
leftphyla<-  c(sort(top10[grep("proteo",top9,invert = T)]),"Others")

# iteration obj.
top10<- c(proteo,leftphyla)
ROW_ORDER<-c()

for(i in top10){
  subdata<- asv_module_mean[asv_module_mean$Phylum==i,1:5]
  p0<- pheatmap(log2(subdata+1),
                cluster_rows=T)
  row_order<- rownames(subdata)[p0$tree_row$order];
  ROW_ORDER<- c(ROW_ORDER,row_order)
}

# got ASV order:
row_annotation<- data.frame(asv_module_mean[ROW_ORDER,"Phylum"])
rownames(row_annotation)<- ROW_ORDER;colnames(row_annotation)="Phylum"

head(row_annotation)

# add colors:
anno_phylum_color<- c("#a04000","#dc7633","#edbb99",
                      "#036EB8","#2A878F","#A48B78",
                      "#b67ccf","#f0d879","#9c9efe","#4F9153")
annorow_color<- list(Phylum=c(Alphaproteobacteria=anno_phylum_color[1],
                              Deltaproteobacteria=anno_phylum_color[2],
                              Gammaproteobacteria=anno_phylum_color[3],
                              Acidobacteria=anno_phylum_color[4],
                              Actinobacteria=anno_phylum_color[5],
                              Bacteroidetes=anno_phylum_color[6],
                              Firmicutes=anno_phylum_color[7],
                              Nitrospirae=anno_phylum_color[8],
                              Verrucomicrobia=anno_phylum_color[9],
                              Others=anno_phylum_color[10]))
pheatmap_gap<- row_annotation;pheatmap_gap$num<-rep(1,nrow(pheatmap_gap))
pheatmap_gap<- pheatmap_gap %>% group_by(Phylum) %>% summarise_each(sum)  %>% data.frame()
rownames(pheatmap_gap)<- pheatmap_gap$Phylum;pheatmap_gap<- pheatmap_gap[top10,-1]
gaps<- c(pheatmap_gap[1])
p3<-pheatmap(log2(asv_module_mean[ROW_ORDER,1:5]+1),
             annotation_row = row_annotation,
             #annotation_col = annotate_col,
             annotation_colors = annorow_color,
             cluster_rows=F,cluster_cols = T,
             color=c(colorRampPalette(colors = viridis(10)[c(1,5)])(10), 
                     colorRampPalette(colors = viridis(10)[c(5:6)])(10),
                     colorRampPalette(colors = viridis(10)[c(6,10)])(100)),
             legend_breaks = c(0,2,4,6,8,10,12),legend_labels=2^c(0,2,4,6,8,10,12))

write.csv(asv_module_mean, paste0(output.path,"Fig2A.heatmap_asv_module_mean.csv"))
save_pheatmap_pdf(p3,paste0(output.path,"Fig2A.Heatmap_asv_5307_5part_mean_log2+1_col2_ann_v2.pdf",
                  width = 4.5,height = 7))

save.image("./1_16S_MH63_NSM_inoculation.RData")
