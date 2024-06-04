#### Script2 #####
# Obj: 16S of MH63_SynCom11, including:
# > Fig2C. Barplot

# URL of Silva138 db
#https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1
#https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1



#### Part I. dada2 pipeline: ASV construction ##########
library(dada2)
library(ggplot2)
library(phyloseq)
library(reshape2)

path<- "~/Documents/Research/LabMember/xjw/Publish/DataAvailibility/Code/"
setwd(path)
output.path<- path
# dada2 pipeline
load("source/2_16S_MH63_SynCom11/2_16S_MH63_dada2.RData")
#rm(errF,errR,mergers,seqtab,seqtab.nochim,taxa,dadaFs,dadaRs,filtFs,filtRs,fnFs,fnRs)

#### Part II. Data preparation ############
library(Biostrings)
library(dplyr)
library(phyloseq)
library(dada2)
library(ggplot2)
library(phyloseq)
library(reshape2)

asv.seq<-colnames(seqtab.nochim)
asv.seq<-DNAStringSet(asv.seq)
names(asv.seq)<-paste0("ASV_",1:length(asv.seq))
#dir.create("2_Process/0_dada2/Table")
#writeXStringSet(asv.seq,"2_Process/0_dada2/Table/Table0.asvseq.fasta")

seqtab.nochim.bk<-seqtab.nochim
colnames(seqtab.nochim.bk)<- paste0("ASV_",1:length(asv.seq))
asv.tab<-t(seqtab.nochim.bk)
rm(seqtab.nochim.bk)
write.csv(asv.tab,"source/2_16S_MH63_SynCom11/Table0.asvtab.csv")

asv.tax<-taxa
rownames(asv.tax)<-paste0("ASV_",1:nrow(asv.tax))
write.csv(asv.tax,"source/2_16S_MH63_SynCom11/Table0.asvtax.csv")


asv.tax<- read.csv("source/2_16S_MH63_SynCom11/Table0.asvtax.csv",header = T,row.names = 1)
asv.tab<- read.csv("source/2_16S_MH63_SynCom11/Table0.asvtab.csv",header = T,row.names = 1)
metadata<- read.table("source/2_16S_MH63_SynCom11/metadata.tsv",sep="\t",header = T)
#metadata<-read.csv("otutab/metadata.csv",header = T,row.names = 1)
rownames(metadata)<- metadata$sampleid

raw.ps<-phyloseq(otu_table(as.matrix(asv.tab),taxa_are_rows = T),
                 tax_table(as.matrix(asv.tax)),
                 sample_data(metadata))

sample_sums(raw.ps)
summary(sample_sums(raw.ps))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 109965  118058  131901  133727  147920  161047 

asv_tax_tab<- cbind(asv.tab, asv.tax)

####### filter
library(dplyr)
library(magrittr)
library(ggplot2)

ps.filt<- prune_taxa(rowSums(otu_table(raw.ps)>0)>=2,raw.ps)
ps.filt<- prune_taxa(data.frame(tax_table(ps.filt))$Kingdom=="Bacteria",ps.filt)
ps.filt<- prune_taxa(data.frame(tax_table(ps.filt))$Family!="Mitochondria",ps.filt) #2943 taxa and 42 samples
colSums(otu_table(ps.filt))[order(colSums(otu_table(ps.filt)))]

filt.ratio<- sample_sums(ps.filt) / sample_sums(raw.ps) %>% as.data.frame()
colnames(filt.ratio)<- "filt.ratio"
filt.ratio$sampleid <- rownames(filt.ratio)
filt.ratio$filtered.depth<- sample_sums(ps.filt) 
filt.ratio$filt.ratio %>% range #0.9536171 0.9992196
p1<- metadata %>%
  inner_join(filt.ratio) %>%
  ggplot(aes(x = sampleid, y = filt.ratio, fill = compartment)) +
  geom_col(width = 0.75) +
  labs(x = "Sample ID", y = "Filtered Ratio", 
       title = "dada2_filt")+
  scale_fill_manual(values = c("#D4C9A8","#6F7179"))+theme_light()+
  theme(axis.text.x = element_text(size=12,color="black",angle=90,vjust=1),
        axis.text.y = element_text(size=12,color="black"),
        axis.title = element_text(size=15, color="black"))+
  labs(y="RA (filtered)", x="",title = "dada2")+
  scale_y_continuous(labels = scales::percent)
p2<- metadata %>%
  inner_join(filt.ratio) %>%
  ggplot(aes(x = sampleid, y = filtered.depth/10000, fill = compartment)) +
  geom_col(width = 0.75) +
  labs(x = "Sample ID", y = "Filtered Ratio", 
       title = "dada2_filt")+
  scale_fill_manual(values = c("#D4C9A8","#6F7179"))+theme_light()+
  theme(axis.text.x = element_text(size=12,color="black",angle=90,vjust=1),
        axis.text.y = element_text(size=12,color="black"),
        axis.title = element_text(size=15, color="black"))+
  labs(y="Depth (filtered) / 10000", x="",title = "dada2")
#ggsave("2_Process/Figure/Figure2.dada2_filtered_Depth.pdf",width = 8,height = 5)

####### check taxonomy
sanger_tax<- readxl::read_excel("source/2_16S_MH63_SynCom11/web.arb-silva.de_align_resultlist_1320026.xlsx",sheet = 2, skip = 1)
colnames(sanger_tax)[1]<- "Stock"
sanger_tax %<>% as.data.frame()
rownames(sanger_tax) <- sanger_tax$Stock

filt_tax<- tax_table(ps.filt)  %>% data.frame()
filt_tab<- otu_table(ps.filt) %>% data.frame()
filt_tab_tax<- cbind(filt_tab, filt_tax)
rm(filt_tax, filt_tab)

## check numbfilt_tab## check number
taxlevel<- c("Phylum","Class","Order","Family","Genus")
tax_stat<- data.frame()
#ht.list<- list()
for(i in 1: length(taxlevel)){
  ta<- taxlevel[i]
  count.sanger <- count(sanger_tax, get(ta));colnames(count.sanger) <- c("Taxon","Ref_sanger")
  count.dada2 <- count(filt_tab_tax, get(ta)); colnames(count.dada2) <- c("Taxon","dada2")
  df<- full_join(count.sanger,count.dada2) %>%replace(is.na(.),0)
  rownames(df)<- df$Taxon
  df$Level<- ta
  tax_stat<- rbind(tax_stat, df)
  # heatmap
  # p<- pheatmap(df[,2:3],show_rownames = T, cluster_rows = F, display_numbers=T, cluster_cols = F, 
  #              fontsize=16,number_format = "%.0f",number_color="black"
  #              #,color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
  #              )
  # ht.list[[i]]<- p
}

## check percentage (RA)
sanger_tax$Family[sanger_tax$Family %in% filt_tab_tax$Family]  
filt_tab_tax$family_present<- ifelse(filt_tab_tax$Family %in% sanger_tax$Family, 1, 0)
#write.csv(filt_tab_tax, "2_Process/0_dada2/dada2.filt_tab_tax.csv")

df<- filt_tab_tax %>% 
  filter(family_present==1) %>% 
  select(starts_with("S11")) %>% 
  t() %>% data.frame() %>%
  mutate(family_reads= rowSums(.))
family_to_filt<- df$family_reads / sample_sums(ps.filt)
family_to_raw <- df$family_reads / sample_sums(raw.ps)

filt.ratio$family_to_filt<- 1-family_to_filt
filt.ratio$family_to_raw<- 1-family_to_raw

index<- c("filt.ratio","family_to_raw","family_to_filt")
filt.ratio$filt.ratio<- 1-filt.ratio$filt.ratio
barlist<- list()
for(i2 in 1:length(index)){
  i<- index[i2]
  p<- metadata %>%
    inner_join(filt.ratio) %>%
    ggplot(aes(x = sampleid, y = get(i), fill = compartment)) +
    geom_col(width = 0.75) +
    labs(x = "Sample ID", y = "Filtered Ratio", 
         title = "dada2_filt")+
    scale_fill_manual(values = c("#D4C9A8","#6F7179"))+theme_light()+
    theme(axis.text.x = element_text(size=12,color="black",angle=90,vjust=1),
          axis.text.y = element_text(size=12,color="black"),
          axis.title = element_text(size=15, color="black"))+
    labs(y="RA (filtered)", x="",title = paste("dada2", i))+
    scale_y_continuous(labels = scales::percent)
  barlist[[i2]]<- p
}
library(gridExtra)
c<-grid.arrange(
  grobs = lapply(barlist, ggplotGrob),nrow=1)
c
#ggsave("2_Process/Figure/Figure3.dada2_barlist.pdf",c,width=18,height = 5)

sanger_tax2<- sanger_tax
sanger_tax$Genus2<- ifelse(is.na(sanger_tax2$Genus),"Rhizobiaceae_Genus",sanger_tax2$Genus)
stock_genus<- sanger_tax[,c("Stock","Genus","Phylum","Class")]
stock_genus$Genus2<- ifelse(is.na(stock_genus$Genus),"Rhizobiaceae_Genus",stock_genus$Genus)
dat<- filt_tab_tax %>%
  mutate(Genus2= ifelse(Family == "Rhizobiaceae", "Rhizobiaceae_Genus",Genus)) %>% 
  mutate(Genus2= ifelse(Genus2 %in% sanger_tax$Genus2, Genus2, "Nontarget")) %>%
  mutate(Genus2= ifelse(is.na(Genus2),"Nontarget",Genus2)) %>% 
  select(starts_with("S11"), Genus2) %>%
  group_by(Genus2) %>% 
  summarise_each(funs = sum) %>% 
  full_join(stock_genus, by="Genus2") %>% replace(is.na(.),"Nontarget")  %>%
  arrange(Phylum,Class,Genus2)
head(dat)
dat[5,"Stock"]<- "W106_ZH1"
dat<- dat[-6,]

dat<- data.frame(dat)
rownames(dat)<- dat$Genus2
a<-dat %>% select(starts_with("S11")) %>%
  apply(2, function(x) x/sum(x)) %>% data.frame() 
dat2<- cbind(a, dat[,c("Stock","Phylum","Class","Genus2")])
dat2$Annotation=paste(dat2$Stock, dat2$Phylum, dat2$Genus2, sep="_") 
dat2$Annotation<- factor(dat2$Annotation, levels=c(dat2$Annotation[-7],dat2$Annotation[7]))


color<- c("#F4A460","#FFCF7F","#F78F8F","#F2D7D5","#8A2BE2","#9C8FBC","#A7D8F2","#7FC7E2","#00BCD4","#A8DADC","#003366")
dat2 %>% 
  reshape2::melt(id.vars=c("Stock","Genus2","Phylum","Class","Annotation")) %>%  
  mutate(Phylum=ifelse(Phylum=="Proteobacteria",Class, Phylum)) %>% 
ggplot(aes(x=variable, y= value, fill=Annotation)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = color)+theme_light()+
  theme(axis.text.x = element_text(size=12,color="black",angle=90,vjust=1),
        axis.text.y = element_text(size=12,color="black"),
        axis.title = element_text(size=15, color="black"))+
  labs(y="RA", x="",title = "dada2")+
  scale_y_continuous(labels = scales::percent)
ggsave("2_Process/Figure/Figure2.dada2_RA.pdf",width = 8,height = 5)


strainlevel<- c("ZH1","W106","ZH398", "W237", "W86", "ZH149", "X137", "W169", "ZH126", "ZH243", "W114")
strainlevel2<- c("W106_ZH1_Firmicutes_Bacillus","ZH398_Firmicutes_Paenibacillus", "W237_Actinobacteriota_Microbacterium", 
                 "W86_Actinobacteriota_Pseudarthrobacter", "ZH149_Proteobacteria_Uliginosibacterium", "X137_Proteobacteria_Burkholderia-Caballeronia-Paraburkholderia", 
                 "W169_Proteobacteria_Pseudomonas", "ZH126_Proteobacteria_Rhizobiaceae_Genus", "ZH243_Bacteroidota_Chitinophaga", "W114_Bacteroidota_Chryseobacterium","Nontarget_Nontarget_Nontarget")
strainlevel3<- c("W169_Proteobacteria_Pseudomonas","W106_ZH1_Firmicutes_Bacillus","W86_Actinobacteriota_Pseudarthrobacter","W237_Actinobacteriota_Microbacterium", "ZH243_Bacteroidota_Chitinophaga",
                 "ZH149_Proteobacteria_Uliginosibacterium","X137_Proteobacteria_Burkholderia-Caballeronia-Paraburkholderia",  "ZH398_Firmicutes_Paenibacillus", "ZH126_Proteobacteria_Rhizobiaceae_Genus",
                 "W114_Bacteroidota_Chryseobacterium")
col<- c("#AC9CC7","#007C87","#F28D8E","#F28D8E","#322C5C","#AC9CC7","#AC9CC7","#007C87","#AC9CC7","#322C5C","#007C87")
dat3<- dat2 %>% 
  reshape2::melt(id.vars=c("Stock","Genus2","Phylum","Class","Annotation")) %>%  
  mutate(Phylum=ifelse(Phylum=="Proteobacteria",Class, Phylum)) %>%  
  mutate(group=ifelse(grepl("S11R",variable),"Root","Soil"))
dat3$Annotation<- factor(dat3$Annotation, levels = strainlevel2 )  
ggplot(dat3, aes(x= Annotation, y=-log10(value), fill=Annotation)) +geom_boxplot(width=.5)+geom_point()+
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1))+ facet_grid( ~ group)+
  theme_light()

#### Part III. plot Figure2C ############
dat4=dat3
dat4$Annotation<- factor(dat4$Annotation, levels = strainlevel3 )  
ggplot(dat4, aes(x= Annotation, y=log10(value), fill=Annotation)) +
  geom_bar(width=.5, stat = "identity")+geom_point(size=.5,alpha=.5)+
  theme_light()+labs(y="log10(RA)")+
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1))+ facet_grid( ~ group)+scale_fill_manual(values = col)
write.csv(dat3,"Fig2C.SynCom11_RA_plotbar.csv")

save.image("2_16S_MH63_SynCom11.RData")
