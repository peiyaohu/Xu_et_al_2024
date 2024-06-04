#### Script3 #####
# Obj: 16S of PRR mutant inoculated with NSM, including:
# > Fig4d-f. diversity
# > Fig4g-i. LEFse visualization: differential taxa (family/genus level)

library(dplyr)
library(tidyverse)
library(Biostrings)
library(gridExtra)
library(agricolae)

#load("/Users/hupeiyao/Documents/Research/LabMember/xjw/16S/xjw_immunityMut/data.RData")
path<- "~/Documents/Research/LabMember/xjw/Publish/DataAvailibility/Code/"
setwd(path)
load("source/3_16S_PRRmut_NSM/3_16S_PRRmut_NSM_dada2.RData")


#### 0. Load Data ############
#asv.seq<-colnames(seqtab.nochim)
#asv.seq<-DNAStringSet(asv.seq)
#names(asv.seq)<-paste0("ASV_",1:length(asv.seq))
#writeXStringSet(asv.seq,"dada2/Table/Table0.asvseq.fasta")

#seqtab.nochim.bk<-seqtab.nochim
#colnames(seqtab.nochim.bk)<- paste0("ASV_",1:length(asv.seq))
#asv.tab<-t(seqtab.nochim.bk)
#rm(seqtab.nochim.bk)
#write.csv(asv.tab,"dada2/Table/Table0.asvtab.csv")

#asv.tax<-readRDS("dada2/RDS/taxa.RData")
#rownames(asv.tax)<-paste0("ASV_",1:nrow(asv.tax))
#write.csv(asv.tax,"dada2/Table/Table0.asvtax.csv")

# asv.tax<- read.csv("dada2/Table/Table0.asvtax.csv",header = T,row.names = 1)
# asv.tab<- read.csv("dada2/Table/Table0.asvtab.csv",header = T,row.names = 1)
# metadata<-read.csv("otutab/metadata.csv",header = T,row.names = 1)
asv.tax<- read.csv("source/3_16S_PRRmut_NSM/Table0.asvtax.csv",header = T,row.names = 1)
asv.tab<- read.csv("source/3_16S_PRRmut_NSM/Table0.asvtab.csv",header = T,row.names = 1)
metadata<-read.csv("source/3_16S_PRRmut_NSM/metadata.csv",header = T,row.names = 1)

raw.ps<-phyloseq(otu_table(as.matrix(asv.tab),taxa_are_rows = T),
                 tax_table(as.matrix(asv.tax)),
                 sample_data(metadata))

summary(colSums(otu_table(raw.ps)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18681   60982   80895   77531   89797  122771
colSums(otu_table(raw.ps))[order(colSums(otu_table(raw.ps)))]

####### filter
ps.filt<- prune_taxa(rowSums(otu_table(raw.ps)>0)>=2,raw.ps)
ps.filt<- prune_taxa(data.frame(tax_table(ps.filt))$Kingdom=="Bacteria",ps.filt)
ps.filt<- prune_taxa(data.frame(tax_table(ps.filt))$Family!="Mitochondria",ps.filt) #2943 taxa and 42 samples
colSums(otu_table(ps.filt))[order(colSums(otu_table(ps.filt)))]

set.seed(26)

ps<-rarefy_even_depth(ps.filt,sample.size = 33742)
ps.removeinput<-prune_samples(data.frame(sample_data(ps))$Compartment!="input",ps) 
ps.removeinput<- prune_taxa(rowSums(otu_table(ps.removeinput))>0,ps.removeinput)  #1917

data<- data.frame(cbind(otu_table(ps.removeinput)))
colnames(data)
geno<- data.frame(Geno=metadata[colnames(data),"Geno"]);rownames(geno)<- colnames(data)
comp<- data.frame(Compartment=metadata[colnames(data),"Compartment"]);rownames(comp)<- colnames(data)
subclass<- rbind(t(geno),t(comp))
write.csv(data,"Fig4d.Sup_Table1.rarefy_33742_tab.csv")
write.csv(tax_table(ps.removeinput),"Fig4d.Sup_Table1.rarefy_33742_tax.csv")

taxonomy<- read.csv("Fig4d.Sup_Table1.rarefy_33742_tax.csv",header = T,row.names = 1,sep=",")

#### 1. Data preparation for LEFse and Visualization (Fig4g-i)####
data2<- as.matrix(data)
data2<- cbind(taxonomy,data2)
data2<- data2 %>% group_by(Taxonomy) %>% 
  summarise_each(funs=sum) %>% data.frame()
rownames(data2)<- data2$Taxonomy;data2<- data2[,-1]
data2<- rbind(subclass,data2)

write.table(data2,"Fig4g_to_h.rarefy_33742_to_LEfSe.txt",sep="\t")

### Visualization: differential taxa (family/genus level)
### boxplot for visualization
boxdata<- read.csv("source/3_16S_PRRmut_NSM/Table3.lefse_all_difftaxa_boxplot.csv",header = T)
boxdata<- boxdata %>%
  filter(boxdata$X =="Keep") 

group<- unique(boxdata$group)
boxlist<- list()
fills<- c("#414141","#18739a","#f7d16c")
for(i in 1:3){
  cg<- group[i]
  subfamily<- boxdata %>%
    dplyr::filter(boxdata$group==cg)
  sub2<- subfamily %>% group_by(variable) %>% summarise_each(mean)
  
  flevel<-  sub2$variable[order(sub2$value,decreasing=T)]
  subfamily$variable<- factor(subfamily$variable,levels = flevel)
  p<- ggplot(subfamily,aes(x = variable, y = value, fill = geno, color = geno)) +
    geom_boxplot(alpha = 0.75, width = .5,outlier.size = .45,size = .25) +
    geom_point(position = position_dodge(width = 0.75), width = 0.5,size=.45) +  # 调整position和width的值
    theme_light() +
    theme(axis.text.x = element_text(size = 12, color = "black",angle = 45,hjust=1),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, color = "black")) +
    xlab("") +ylab("RA")+
    scale_color_manual(values = fills) +
    scale_fill_manual(values = fills) +
    ggtitle(cg)
  boxlist[[i]]<- p
}

c<-grid.arrange(
  grobs = lapply(boxlist, ggplotGrob))
ggsave("Fig4g_h.Figure.lefse_box_keep.pdf",c,width=4,height = 10)


#### 2. Diversity (Fig4d-f) ############
### 1. basis stats
stats<- data.frame(sampleid=rownames(sample_data(raw.ps)))
stats$rawreads<- colSums(otu_table(raw.ps))[stats$sampleid]
stats$bactreads<- colSums(otu_table(ps.filt))[stats$sampleid]
stats$bactratio<- round(stats$bactreads/stats$rawreads*100,2)
rownames(stats)<-stats$sampleid
View(stats)

### 2. alpha diversity (Fig4d)
sample_sums(ps)
alpha<- estimate_richness(ps,measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))

metadata_alpha<- cbind(metadata, alpha[rownames(metadata),])
stats<- cbind(stats, metadata_alpha[rownames(stats),])
write.csv(stats,"Fig4d_Sup_NGS_stats_metadata_alpha.csv")

color7<- c("#263859","#3b8686","#79bd9a","#f06966","#f1ac9d","#1F6ED4","#39BAE8")
color7_2<-color7[c(1,2,4,6,3,5,7)]
method = c("Observed","Shannon","Simpson","Chao1")
sub_design<- metadata_alpha
level1<- c("input","Ksoil","rsoil","csoil","kroot","rroot","croot")
sub_design$Group<- factor(sub_design$Group, levels = level1)
sub_design<- sub_design %>%
  filter(rownames(sub_design) %in% sample_names(ps))


for(m in method){
  #m = "shannon"
  model = aov(sub_design[[m]] ~ Group, data=sub_design)
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  Tukey_HSD_table = as.data.frame(Tukey_HSD$Group) 
  write.table(paste(m, "\n\t", sep=""), file=paste("Fig4_Tab_alpha_",m,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  suppressWarnings(write.table(Tukey_HSD_table, file=paste("Fig4_Tab_alpha_",m,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
  
  # LSD test for stat label
  out = LSD.test(model,"Group", p.adj="fdr") # alternative fdr
  stat = out$groups
  sub_design$stat=stat[as.character(sub_design$Group),]$groups
  max=max(sub_design[,c(m)])
  min=min(sub_design[,c(m)])
  x = sub_design[,c("Group",m)]
  y = x %>% group_by(Group) %>% summarise_(Max=paste('max(',m,')',sep=""))
  y=as.data.frame(y)
  rownames(y)=y$Group
  sub_design$y=y[as.character(sub_design$Group),]$Max + (max-min)*0.1
  p<-ggplot(sub_design,aes(x=Group, y=sub_design[[m]]))+
    geom_boxplot(aes(fill=Group),outlier.size = .5,linewidth=.3)+
    #scale_fill_manual(
    #values=c("#00BFC4", "#F9766E", "#9E7C61", "darkorange"))+
    geom_jitter(position=position_jitter(0.17), size=.5, alpha=0.7)+scale_y_continuous(limits = c(0,1.2*max))+
    labs(title="",x="Compartment",y=m)+
    geom_text(data=sub_design, aes(x=Group, y=y, color=Group, label= stat))+
    theme(
      legend.position="none",
      axis.title.y=element_text(hjust=.5),
      plot.title=element_text(size=11,hjust=.5),
    )+theme_bw()+scale_fill_manual(values = color7_2)+scale_color_manual(values = color7_2)+
    theme(
      axis.text.x = element_text(size=12,angle = 45,hjust = 1,color="black"),
      axis.text.y = element_text(size=12,color="black"),
      axis.title = element_text(size = 15,color="black")
    )
  ggsave(paste("Fig4_", m, ".pdf", sep=""), p, width = 10, height = 9,units = "cm")
  p
}         

#### 3. CPCoA (Fig4e-f) ##### 
library(amplicon)
# set factor for sample_data
sample_data_ps<- data.frame(sample_data(ps)) %>%
  rownames_to_column(var="ID")
sample_data_ps$Geno=factor(sample_data_ps$Geno, levels = c("input","Kitaake","fls2","cerk1"))

asvtab_r<- ps %>% otu_table() %>% data.frame() 
beta<- vegdist(t(asvtab_r),method = "bray")

beta_df<- beta %>%  
    cmdscale(k=4) %>% as.data.frame()%>% rownames_to_column(var="ID")
eig<- beta %>%  
    cmdscale(k=4, eig=T) %>% .$eig

colnames(beta_df)<- c("ID",paste0("PC",1:4))
beta_df<- inner_join(beta_df,sample_data_ps) %>% as_tibble()
beta_df$Group<- factor(beta_df$Group, levels=c("input", "Ksoil","rsoil","csoil",
                                                 "kroot","rroot","croot"))
centroid<-beta_df %>%
    select(starts_with("PC"),Group) %>%
    group_by(Group) %>%
    summarise(PC1=mean(PC1),
              PC2=mean(PC2))
# PCoA
ggplot(beta_df, aes(x=PC1, y=PC2, color=Group))+
    geom_point(alpha=.7, size=2) +
    labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
         y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
         title=paste(m," PCoA",sep=""))+
    geom_point(data=centroid,aes(fill=Group),color="black",shape=21,size=2,show.legend = F)+
    scale_fill_manual(values = color7_2)+
    scale_color_manual(values = color7_2)+
    theme_light()+
    theme(axis.title = element_text(size=15,colour = "black"),
          axis.text = element_text(size=12,colour = "black"),
          legend.text= element_text(size=6,colour = "black"),
          legend.title = element_text(size=7),
          legend.key.size = unit(.3, 'cm'),
          legend.margin = margin(t=1,r=1,b=1,l=1,unit="pt"))+stat_ellipse(show.legend = F)

# PERMANOVA
perm<-1e4
permanova_df<- data.frame(method=c(),
                          obj=c(),
                          Df=c(),
                          SumOfSqs=c(),
                          R2=c(),
                          F=c(),
                          p_value=c(),
                          permutations=c())
permanova_non_df<- data.frame(method=c(),
                              obj=c(),
                              Df=c(),
                              SumOfSqs=c(),
                              R2=c(),
                              F=c(),
                              p_value=c(),
                              permutations=c())

main_theme =   theme(axis.line.x=element_line(size=.5, colour="black"),
                     axis.line.y=element_line(size=.5, colour="black"),
                     axis.ticks=element_line(color="black"),
                     axis.text=element_text(color="black", size=12))
permanova_non<- adonis2(beta~beta_df$Compartment+beta_df$Geno,permutations = perm)
stats_df_non<- data.frame(method=rep(m,4),
                            obj=c("Compartment","Geno","Residual","Total"),
                            Df=permanova_non$Df,
                            SumOfSqs=permanova_non$SumOfSqs,
                            R2=permanova_non$R2,
                            p_value=permanova_non$`Pr(>F)`,
                            permutations=rep(perm,4)) 
permanova_non_df<- rbind(permanova_non_df,stats_df_non)
  
permanova<- adonis2(beta~beta_df$Compartment*beta_df$Geno,permutations = perm)
stats_df<- data.frame(method=rep(m,5),
                        obj=c("Compartment","Geno","Interaction","Residual","Total"),
                        Df=permanova$Df,
                        SumOfSqs=permanova$SumOfSqs,
                        R2=permanova$R2,
                        F=permanova$F,
                        p_value=permanova$`Pr(>F)`,
                        permutations=rep(perm,5)) 
permanova_df<- rbind(permanova_df,stats_df)
  
# Compartment dispersions
comparison_groups<- list(c("Ksoil","rsoil","csoil"),
                           c("kroot","rroot","croot"))
names(comparison_groups)<- c("soil","root")
for(i in 1:2){
    sub_sample<- sample_data_ps %>%
      filter(Group %in% comparison_groups[[i]]) 
    sub_beta<- beta %>% as.matrix() %>%
      .[sub_sample$ID,sub_sample$ID] %>%
      as.dist()
    
    sub_sample$Geno<- factor(sub_sample$Geno,levels = c("Kitaake","fls2","cerk1"))
    sub_tab<- asvtab_r[,sub_sample$ID] 
    rownames(sub_sample)<- sub_sample$ID
    if(m!="aitchison"){
      p<-beta_cpcoa(sub_tab,sub_sample,dis=m,groupID = "Geno",ellipse = T,label=F)+
        scale_color_manual(values = c("#18739a","#414141","#f7d16c"))+main_theme
      p
      ggsave(paste0("Fig4e_to_f.CPCoA_",names(comparison_groups)[i],"_",m,".pdf"),width = 4,height = 3)
    }
}  
    
save.image("3_16S_PRRmut_NSM.RData")






