#The script is for Figure 2
# Date: 2024/07/05
# By: hupeiyao.


#### Part I. ASV construction ##########
library(dada2)
library(ggplot2)
library(phyloseq)
library(reshape2)

path<- "~/Documents/Research/LabMember/xjw/16S/20240719/"
setwd(path)
rawseq.dir<-paste0(path,"1_cutadapt")
fnFs <- sort(list.files(rawseq.dir, pattern="-R1.fq.gz.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(rawseq.dir, pattern="-R2.fq.gz.trimmed.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)
QC.F<-plotQualityProfile(fnFs[1:6])
QC.R<-plotQualityProfile(fnRs[1:6])
ggsave("2_Process/Figure/QC_F.pdf",QC.F,width=20,height=20,units = "cm")
ggsave("2_Process/Figure/QC_R.pdf",QC.R,width=20,height=20,units = "cm")

filtFs <- file.path(path, "2_Process/0_dada2/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "2_Process/0_dada2/filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)# On Windows set multithread=FALSE
head(out)
write.csv(out,"2_Process/0_dada2/filt_statout.csv")

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
saveRDS(mergers,"2_Process/0_dada2/mergers.RData")
saveRDS(dadaFs,"2_Process/0_dada2/dadaFs.RData")
saveRDS(dadaRs,"2_Process/0_dada2/dadaRs.RData")
saveRDS(errF,"2_Process/0_dada2/errF.RData")
saveRDS(errR,"2_Process/0_dada2/errR.RData")
saveRDS(filtFs,"2_Process/0_dada2/filtFs.RData")
saveRDS(filtRs,"2_Process/0_dada2/filtRs.RData")

seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)#12 6835

sum(seqtab.nochim)/sum(seqtab) #0.7949234

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, "2_Process/0_dada2/dada2_track.csv")

head(track)

seqtab_file <- "~/Documents/Research/github_16S/silva138/silva_nr99_v138.1_train_set.fa"
species_file <- "~/Documents/Research/github_16S/silva138/silva_species_assignment_v138.1.fa"
# 
# batch_size <- 1000
# 
# num_batches <- ceiling(ncol(seqtab.nochim) / batch_size)
# 
# taxa_list <- list()
# 
# for (i in 1:num_batches) {
#   start_idx <- (i - 1) * batch_size + 1
#   end_idx <- min(i * batch_size, nrow(seqtab.nochim))
#   
#   seqtab_subset <- seqtab.nochim[start_idx:end_idx, ]
#   
#   taxa_subset <- assignTaxonomy(seqtab_subset, seqtab_file, multithread=TRUE)
#   taxa_subset <- addSpecies(taxa_subset, species_file)
#   
#   taxa_list[[i]] <- taxa_subset
# }
# 
# taxa_combined <- do.call(rbind, taxa_list)
# saveRDS(taxa_combined, "taxa.rds")
# dim(taxa_combined)

saveRDS(seqtab.nochim, "2_Process/0_dada2/seqtab.nochim.RData")

#rm(errF,errR,mergers,seqtab,seqtab.nochim,taxa,dadaFs,dadaRs,filtFs,filtRs,fnFs,fnRs)
save.image(file = "2_Process/0_dada2/dada2.RData")
#### Part II. Phyloseq-object construct ############
library(Biostrings)
library(dplyr)
library(phyloseq)
library(dada2)
library(ggplot2)
library(phyloseq)
library(reshape2)

path<- "~/Documents/Research/LabMember/xjw/16S/20240719/"
setwd(path)
# load("2_Process/0_dada2/dada2.RData")
# 
# asv.seq<-colnames(seqtab.nochim)
# asv.seq<-DNAStringSet(asv.seq)
# names(asv.seq)<-paste0("ASV_",1:length(asv.seq))
#dir.create("2_Process/0_dada2/Table")
# writeXStringSet(asv.seq,"2_Process/Table/Table0.asvseq.fasta")

# seqtab.nochim.bk<-seqtab.nochim
# colnames(seqtab.nochim.bk)<- paste0("ASV_",1:length(asv.seq))
# asv.tab<-t(seqtab.nochim.bk)
# rm(seqtab.nochim.bk)
# write.csv(asv.tab,"2_Process/Table/Table0.asvtab.csv")
# 
# #asv.tax<-readRDS("2_Process/0_dada2/taxa2.RData")
# asv.tax<- taxa
# rownames(asv.tax)<-paste0("ASV_",1:nrow(asv.tax))
# write.csv(asv.tax,"2_Process/Table/Table0.asvtax.csv")

asv.tax<- read.csv("2_Process/Table/Table0.asvtax.csv",header = T,row.names = 1)
asv.tab<- read.csv("2_Process/Table/Table0.asvtab.csv",header = T,row.names = 1)
metadata<- read.table("source/1.sample-metadata.csv",sep=",",header = T)
#metadata<-read.csv("otutab/metadata.csv",header = T,row.names = 1)
rownames(metadata)<- metadata$sampleid

raw.ps<-phyloseq(otu_table(as.matrix(asv.tab),taxa_are_rows = T),
                 tax_table(as.matrix(asv.tax)),
                 sample_data(metadata))

sample_sums(raw.ps)
summary(sample_sums(raw.ps))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9897   57742   75985   72892   88150  140194 

asv_tax_tab<- cbind(asv.tab, asv.tax)

save.image("2_Process/0_dada2/raw.ps.RData")

#### Part III. Diversity (Figure 2A-C)####
# Subset WT (MH63 and Kitaake)
sampledata<- data.frame(sample_data(raw.ps))
two_WT_ps<- prune_samples(sampledata$Geno %in% c("Kitaake", "MH63", "Input") & sampledata$SoilType == "Loam", raw.ps) 
rm(sampledata)

source("~/Documents/Research/mGWAS/mGWAS2024/999_R_function/ps_filter.R")
sample_data(two_WT_ps) %>% data.frame() %>% count(Batch, Compartment, Geno)

# Basic control: remove taxa occurred less than 6 samples.
ps.filter<- ps_filter(two_WT_ps, prev_n = 6, outdir = "./2_Process/Table/")
ps.filt<- ps.filter$ps.filter #1104 taxa and 36 samples 
write.csv(tax_table(ps.filt), "2_Process/Table/Table0.asvtax_WT_associated.csv")
write.csv(otu_table(ps.filt), "2_Process/Table/Table0.asvtab_WT_associated.csv")

filt.ratio<- ps.filter$discard_log
sampledata<- data.frame(sample_data(ps.filt))
filt.ratio %>% 
  select(!percentage) %>%
  mutate(low_prev_pct= LowPrevalence/Raw,
         kingdom_pct= Kingdom/Raw,
         mitochondria_pct = Mitochondria/Raw,
         Retained_pct = Retained/Raw) %>%
  select(ends_with("pct")) %>% 
  tibble::rownames_to_column("sampleid") %>% 
  left_join(sampledata[,c("Geno","sampleid")],by="sampleid") %>% 
  select(!sampleid) %>% 
  melt(id.vars = "Geno") %>% 
  group_by(variable, Geno) %>% 
  summarise(mean_value = mean(value, na.rm = TRUE)) %>% 
  ggplot(aes(x=Geno, y=mean_value, fill=variable)) +
  geom_bar(stat = "identity", position = "stack", width=.60)+
  scale_fill_manual(values=c("#008080","#5DAE8B","#0D98BA","#ADD8E6"))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black",size=12, angle=90, hjust = 1),
        axis.text.y = element_text(color="black",size=12),
        axis.title = element_text(color="black",size=15)) +
  labs(title = "Basic control")+scale_x_discrete(labels = c("All Input (12)", "Kitaake Root-associated (12)", "MH63 Root-associated (12)"))
ggsave("2_Process/Figure/Figure1.filter.pdf", width=6, height = 8)


### Diversity
source("~/Documents/Research/mGWAS/mGWAS2024/999_R_function/alpha_div.R")
sampledata<- data.frame(sample_data(ps.filt))
sampledata %<>% mutate(Batch= ifelse(Batch %in% c(1,2),1,2)) %>%
  mutate(Group= paste(Geno, Compartment, Batch, sep="_")) %>%
  mutate(Group= factor(Group, levels= c("Input_Input_1","MH63_Rhizo_1", "MH63_Root_1", "Input_Input_2", "Kitaake_Rhizo_2", "Kitaake_Root_2"),
                       labels=c("MH63_Input","MH63_Rhizo", "MH63_Root", "Kitaake_Input", "Kitaake_Rhizo", "Kitaake_Root")))
sampledata$Group

colorRGB<- c("#795548","#1A237E","#006064",
             "#D4938A","#3F51B5","#00BCD4")
group_levels = c("MH63_Input","MH63_Rhizo", "MH63_Root", "Kitaake_Input", "Kitaake_Rhizo", "Kitaake_Root")

## Figure2A-B
alpha_div(ps = ps.filt, metadata= sampledata, outdir="2_Process/",
          group_levels = group_levels,
          colorRGB = colorRGB, xlab_title = "Group",  tabidx = 1, figidx = 2)

## Figure2C
source("~/Documents/Research/mGWAS/mGWAS2024/999_R_function/beta_div.R")
beta_div(ps.filt, sample_data_ps=sampledata, colname="sampleid", group_levels= group_levels,
         colorRGB=colorRGB, outdir="2_Process/", tabidx=1, figidx=2, perm_group = c("Compartment","Geno","Batch"))


#### Part IV. core microbiota (Figure 2D-E) ####

retain_prevalence <- function(ps, tabind = 1, figind = 1, filtN = 0, RA = 0, prevN = 0, 
                              var_name = NULL, out_dir = "./", grid_ncol = 1, grid_nrow = 3,grid_width=12, grid_height=12) {
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  
  if (!dir.exists(file.path(out_dir, "Figure"))) { dir.create(file.path(out_dir, "Figure")) }
  if (!dir.exists(file.path(out_dir, "Table"))) { dir.create(file.path(out_dir, "Table")) }
  
  # A function to get otutab from phyloseq obj.
  get_otutab <- function(ps) { data.frame(otu_table(ps)) }
  
  ######## required parameters ######
  # filtN (integer)
  # tabind, figind (integer, output index)
  # RA (0 ~ 1)
  # prevN (0 ~ 1)
  # var_name:
  ##########################
  
  # define the theme
  hupy.theme = theme(axis.text.x = element_text(size = unit(7, "pt")),
                     axis.line.x = element_line(size = .5, colour = "black"),
                     axis.line.y = element_line(size = .5, colour = "black"),
                     axis.ticks = element_line(color = "black"),
                     axis.text = element_text(color = "black", size = unit(7, "pt")),
                     axis.title = element_text(size = unit(7, "pt")),
                     legend.text = element_text(size = 6, color = "black"),
                     legend.title = element_text(size = 7, colour = "black"),
                     legend.key.size = unit(.3, 'cm'),
                     legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"))
  
  process_ps <- function(subps, varname) {
    rawotu <- get_otutab(subps)
    otu_table <- rawotu[rowSums(rawotu) >= filtN,]
    depth <- colSums(rawotu)
    
    # calculate prevalence of each OTU
    prevalence <- rowSums(rawotu > 0) / ncol(rawotu)
    # calculate RA of OTUs
    otu.ra <- as.data.frame(apply(rawotu, 2, function(x) x / sum(x)))
    otu.ra$mean <- apply(otu.ra, 1, mean)
    
    # Combine prevalence and RA together
    prev.ra <- data.frame(ID = rownames(otu.ra),
                          prev = prevalence,
                          RA.mean = otu.ra$mean)
    prev.ra$log.mean <- log10(prev.ra$RA.mean)
    
    # add conditional column (filtN、ra、prevN) and set color.
    prev.ra <- prev.ra %>%
      mutate(
        filtN = ifelse(ID %in% rownames(otu_table), "Retain", "Discard"),
        ra = ifelse(RA.mean > RA, "Retain", "Discard"),
        prevN = ifelse(prev > prevN, "Retain", "Discard"),
        condition = case_when(
          ra == "Retain" & prevN == "Retain" ~ "Both_Pass",
          ra == "Discard" & prevN == "Discard" ~ "Both_Fail",
          ra == "Retain" & prevN == "Discard" ~ "Only_RA_pass",
          ra == "Discard" & prevN == "Retain" ~ "Only_prevN_pass",
          TRUE ~ ""
        ),
        condition = factor(condition, levels = c("Both_Pass", "Only_RA_pass", "Only_prevN_pass", "Both_Fail", ""))
      )
    write.csv(prev.ra, paste(out_dir, "Table/Table", tabind, ".", varname, "_prev", prevN, "_RA", RA, "_filtN", filtN, ".csv", sep = ""))
    
    RA_stats <- prev.ra %>%
      dplyr::select(condition, RA.mean) %>%
      dplyr::group_by(condition) %>%
      dplyr::summarise(RA_mean_sum = sum(RA.mean)) %>%
      left_join(prev.ra %>% count(condition), by = "condition") %>%
      mutate(variable = varname)
    
    pass_otu_ra <- function(data, para_col, value_col) {
      data %>%
        filter({{ para_col }} == "Retain") %>%
        summarise(total = sum({{ value_col }}, na.rm = T)) %>%
        dplyr::pull(total) %>% round(3) # Extract a single column
    }
    
    process <- data.frame(
      parameter = c("raw", "filtN", "RA", "prev"),
      threshold = c("-", filtN, RA, prevN),
      n_pass_otu = c(nrow(prev.ra),
                     length(which(prev.ra$filtN == "Retain")),
                     length(which(prev.ra$ra == "Retain")),
                     length(which(prev.ra$prevN == "Retain"))),
      n_pass_otu_ra = c(1,
                        pass_otu_ra(prev.ra, filtN, RA.mean),
                        pass_otu_ra(prev.ra, ra, RA.mean),
                        pass_otu_ra(prev.ra, prevN, RA.mean)
      ),
      variable = varname
    )
    
    stats_list <- list(RA_stats, process)
    names(stats_list) <- c("RA_stats", "process")
    
    # plot
    p <- ggplot(prev.ra, aes(x = log.mean, y = prev, color = condition, fill = condition)) + 
      geom_point(alpha = .8) +
      theme_classic() + hupy.theme +
      labs(x = "log(Mean_RelativeAbundance)", y = "Prevalence", title = varname) +
      scale_color_manual(values = c("#0b4d4d","#7FB3D5","#A8DADC","#DEDEDE", "#C89770")) +
      geom_vline(xintercept = log10(RA), linetype = "dashed", color = "grey") +
      geom_hline(yintercept = prevN, linetype = "dashed", color = "grey") +
      # annotate("text", x = -6, y = 0.875,
      #          label = paste("RA of retained OTUs:", round(RA_stats$RA_mean_sum[1], 3)),
      #          size = 2.5) +
      # annotate("text", x = -6, y = 0.75,
      #          label = paste("RA = ", RA, "; filtN = ", filtN, "; prevN = ", prevN, sep = " "), size = 2.5) +
      ylim(0, 1) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 15))
    ggsave(paste(out_dir, "Figure/Figure", figind, ".", varname, ".RA_", RA, "_prevN_", prevN, ".pdf", sep = ""), plot = p, 
           width = 4, height = 3 )
    return(list(stats_list = stats_list, plot = p))
  }
  
  if (inherits(ps, "list")) {
    all_stats <- lapply(names(ps), function(varname) process_ps(ps[[varname]], varname))
    
    RA_stats_df_total <- do.call(rbind, lapply(all_stats, function(x) x$stats_list$RA_stats))
    process_df_total <- do.call(rbind, lapply(all_stats, function(x) x$stats_list$process))
    
    plots <- lapply(all_stats, function(x) x$plot)
    title <- ggdraw() + draw_label("OTU retained under different standards", fontface = 'bold')
    final_plot <- plot_grid(plotlist = plots, ncol = grid_ncol, nrow = grid_nrow, labels = "AUTO")
    final_plot_with_title <- plot_grid(title, final_plot, nrow = 2, rel_heights = c(0.2, 1))
    ggsave(paste(out_dir, "Figure/Grid_", var_name, ".pdf", sep = ""), plot = final_plot_with_title, width = grid_width, height = grid_height)
    
    return(list(RA_stats_df_total = RA_stats_df_total, process_df_total = process_df_total))
  } else {
    if (is.null(var_name)) {
      stop("var_name must be provided when ps is not a list")
    }
    result <- process_ps(ps, var_name)
    return(result$stats_list)
  }
}

sample_data(ps.filt) %>% data.frame() %>% count(Batch, Compartment, Geno)


ps_mh_rhizo<- prune_samples(data.frame(sample_data(ps.filt))$Compartment == "Rhizo" & data.frame(sample_data(ps.filt))$Geno == "MH63", ps.filt)
ps_mh_root<- prune_samples(data.frame(sample_data(ps.filt))$Compartment == "Root" & data.frame(sample_data(ps.filt))$Geno == "MH63", ps.filt)
ps_kit_rhizo<- prune_samples(data.frame(sample_data(ps.filt))$Compartment == "Rhizo" & data.frame(sample_data(ps.filt))$Geno == "Kitaake", ps.filt)
ps_kit_root<- prune_samples(data.frame(sample_data(ps.filt))$Compartment == "Root" & data.frame(sample_data(ps.filt))$Geno == "Kitaake", ps.filt)

ps_list_1<- list(ps_mh_rhizo= ps_mh_rhizo,
                 ps_mh_root= ps_mh_root,
                 ps_kit_rhizo= ps_kit_rhizo,
                 ps_kit_root= ps_kit_root)

list_1<- retain_prevalence(ps_list_1, filtN = 0,  RA = 0.01, prevN = 5/6, 
                           out_dir = "~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/", tabind = 1, figind = 1,grid_ncol = 4,
                           grid_width = 18, grid_height = 12,var_name = "ra0.01")

list_1$RA_stats_df_total %>%
  write.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table_conditionsummary1_res_RA0.01_prev0.83.csv")

list_1$process_df_total%>%
  write.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table_conditionsummary2_res_RA0.01_prev0.83.csv")


# res_mh_rhizo<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/Table/Table1.ps_mh_rhizo_prev0.833333333333333_RA0.001_filtN0.csv", header = T, row.names = 1)
# res_mh_root<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/Table/Table1.ps_mh_root_prev0.833333333333333_RA0.001_filtN0.csv", header = T, row.names = 1)
# res_kit_rhizo<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/Table/Table1.ps_kit_rhizo_prev0.833333333333333_RA0.001_filtN0.csv", header = T, row.names = 1)
# res_kit_root<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/Table/Table1.ps_kit_root_prev0.833333333333333_RA0.001_filtN0.csv", header = T, row.names = 1)

res_mh_rhizo<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table1.ps_mh_rhizo_prev0.833333333333333_RA0.01_filtN0.csv", header = T, row.names = 1)
res_mh_root<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table1.ps_mh_root_prev0.833333333333333_RA0.01_filtN0.csv", header = T, row.names = 1)
res_kit_rhizo<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table1.ps_kit_rhizo_prev0.833333333333333_RA0.01_filtN0.csv", header = T, row.names = 1)
res_kit_root<- read.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table1.ps_kit_root_prev0.833333333333333_RA0.01_filtN0.csv", header = T, row.names = 1)

colnames(res_mh_rhizo)<- paste("mh_rhizo",colnames(res_mh_rhizo),sep="_")
colnames(res_mh_root)<- paste("mh_root",colnames(res_mh_root),sep="_")
colnames(res_kit_rhizo)<- paste("kit_rhizo",colnames(res_kit_rhizo),sep="_")
colnames(res_kit_root)<- paste("kit_root",colnames(res_kit_root),sep="_")
colnames(res_mh_rhizo)[1]<- "ASVID"
colnames(res_mh_root)[1]<- "ASVID"
colnames(res_kit_rhizo)[1]<- "ASVID"
colnames(res_kit_root)[1]<- "ASVID"

library(tidyverse)
asv_tax<- ps.filt %>% 
  tax_table() %>%
  data.frame() 
all_res<- res_mh_rhizo %>%
  left_join(res_mh_root) %>% 
  left_join(res_kit_rhizo) %>% 
  left_join(res_kit_root) %>%
  left_join(asv_tax %>% rownames_to_column("ASVID")) 
write.csv(all_res, "~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table_all_res_RA0.01_prev0.83.csv")
all_res %>%
  # select(ends_with("condition"),ASVID,Phylum,Class,Order,Family,Genus) %>%
  filter((mh_rhizo_condition == "Both_Pass"&mh_root_condition == "Both_Pass")&
           (kit_rhizo_condition == "Both_Pass"&kit_root_condition == "Both_Pass") )  %>%
  write.csv("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/Table/Table_filt_res_RA0.01_prev0.83.csv")


# filter
WT_ps<- prune_samples(data.frame(sample_data(ps.filt))$Geno %in% c("MH63","Kitaake"), ps.filt)
WT_ps_ra<- WT_ps %>% otu_table() %>% data.frame() %>% 
  apply(2, function(x) x/sum(x)) %>% data.frame
thred_ra<- 1e-3

WT_ps_ra %>% apply(2, function(x) sum(x>1e-2)) %>% data.frame() %>% write.csv("asv_pass_thred1.csv")
asv_pass_df<- WT_ps_ra %>% apply(2, function(x) x>thred_ra) %>%
  as.data.frame() 
asv_pass<- asv_pass_df%>%
  apply(1, function(x) sum(x==1)) %>% data.frame
colnames(asv_pass)<- "N"
asvid<- asv_pass %>% filter(N>0) %>%rownames_to_column("ASVID") %>% pull(ASVID)
WT_ps_ra[asvid,] %>%
  rownames_to_column("ASVID") %>%
  left_join(asv_pass %>% filter(N>0) %>% rownames_to_column("ASVID")) %>% 
  left_join(data.frame(tax_table(WT_ps)) %>% rownames_to_column("ASVID")) %>% 
  write.csv(paste0("WT_ps_ra_individual_",thred_ra, ".csv"))


###### redraw: highligh ASV_12
hupy.theme = theme(axis.line.x = element_line(size = .5, colour = "black"),
                   axis.line.y = element_line(size = .5, colour = "black"),
                   axis.ticks = element_line(color = "black"),
                   axis.text = element_text(color = "black", size = unit(7, "pt")),
                   axis.title = element_text(size = unit(7, "pt")),
                   legend.text = element_text(size = 6, color = "black"),
                   legend.title = element_text(size = 7, colour = "black"),
                   legend.key.size = unit(.3, 'cm'),
                   legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"))


RA=0.01
prevN = 5/6




######## Venn #
all_res<- res_mh_rhizo %>%
  left_join(res_mh_root) %>% 
  left_join(res_kit_rhizo) %>% 
  left_join(res_kit_root) %>%
  left_join(asv_tax %>% rownames_to_column("ASVID")) 

group_pass_id<- list(MH63_Rhizo = all_res %>%  filter(mh_rhizo_condition == "Both_Pass") %>% pull(ASVID),
                     MH63_Root = all_res %>%  filter(mh_root_condition == "Both_Pass") %>% pull(ASVID),
                     Kit_Rhizo = all_res %>%  filter(kit_rhizo_condition == "Both_Pass") %>% pull(ASVID),
                     Kit_Root = all_res %>%  filter(kit_root_condition == "Both_Pass") %>% pull(ASVID))
ggvenn::ggvenn(group_pass_id)

colorRGB<- c("#795548","#1A237E","#006064",
             "#D4938A","#3F51B5","#00BCD4")

venn<- ggvenn::ggvenn(group_pass_id,
                      fill_color = c(colorRGB[c(2,3,5,6)]))
ggsave("Figure2.Venn.pdf", width=5, height=5)


##### Venn2
all_res<- all_res %>%
  mutate(mh_rhizo_present = ifelse(mh_rhizo_RA.mean == 0, "Absent", "Present"),
         mh_root_present = ifelse(mh_root_RA.mean ==0, "Absent", "Present"),
         kit_rhizo_present = ifelse(kit_rhizo_RA.mean == 0,"Absent", "Present"),
         kit_root_present = ifelse(kit_root_RA.mean ==0, "Absent", "Present"))

group_present_id<- list(MH63_Rhizo = all_res %>%  filter(mh_rhizo_present == "Present") %>% pull(ASVID),
                        MH63_Root = all_res %>%  filter(mh_root_present == "Present") %>% pull(ASVID),
                        Kit_Rhizo = all_res %>%  filter(kit_rhizo_present == "Present") %>% pull(ASVID),
                        Kit_Root = all_res %>%  filter(kit_root_present == "Present") %>% pull(ASVID))

colorRGB<- c("#795548","#1A237E","#006064",
             "#D4938A","#3F51B5","#00BCD4")

venn2<- ggvenn::ggvenn(group_present_id,
                       fill_color = c(colorRGB[c(2,3,5,6)]))
ggsave("Figure2.Venn2.pdf", width=5, height=5)


#### PartV. pNST (Figure 2F) ####

path<- "~/Documents/Research/LabMember/xjw/16S/20240719/"
setwd(path)
load("2_Process/0_dada2/raw.ps.RData")
taxa<- readRDS("~/Documents/Research/LabMember/xjw/16S/20240719/2_Process/0_dada2/taxa.RData")
# asv_seq<- Biostrings::DNAStringSet(rownames(taxa))
# names(asv_seq)<- paste0("ASV_",1:length(asv_seq))
# Biostrings::writeXStringSet(asv_seq, "~/Documents/Research/LabMember/xjw/16S/20240719/2_Process/0_dada2/asv_seq_20250207.fa")
asv_seq<- Biostrings::readDNAStringSet("~/Documents/Research/LabMember/xjw/16S/20240719/2_Process/0_dada2/asv_seq_20250207.fa")

raw.ps #2627AA_ALPHABETraw.ps #26277 taxa taxaraw.ps #26277 taxa and 72 samples

# Subset WT (MH63 and Kitaake)
sampledata<- data.frame(sample_data(raw.ps))
two_WT_ps<- prune_samples(sampledata$Geno %in% c("Kitaake", "MH63", "Input") & sampledata$SoilType == "Loam", raw.ps) 
rm(sampledata)

source("~/Documents/Research/mGWAS/mGWAS2024/999_R_function/ps_filter.R")
sample_data(two_WT_ps) %>% data.frame() %>% count(Batch, Compartment, Geno)

# Basic control: remove taxa occurred less than 6 samples.
ps.filter<- ps_filter(two_WT_ps, prev_n = 6, outdir = "./2_Process/Table/")

####### Stochastic process/Deterministic process
### Date：2025/2/7
asv_seq<- Biostrings::readDNAStringSet("~/Documents/Research/LabMember/xjw/16S/20240719/2_Process/0_dada2/asv_seq_20250207.fa")
ps.filt<- prune_samples(data.frame(sample_data(ps.filter$ps.filter))$Geno %in% c("MH63","Kitaake"), ps.filter$ps.filter)

ps.filt_seq<- asv_seq[taxa_names(ps.filt)]
# writeXStringSet(ps.filt_seq, "20250207/ps.filt_seq_1104.fa")
##### conduct in qiime2 for phylogenitic tree
# /public/home/pyhu/project/xjw/16S/20250207_phylo
# qiime phylogeny align-to-tree-mafft-fasttree \
# --i-sequences ps.filt_seq_1104.qza \
# --o-alignment 1_aligned-rep-seqs.qza \
# --o-masked-alignment 2_masked-aligned-rep-seqs.qza \
# --o-tree 3_unrooted-tree.qza \
# --o-rooted-tree 4_rooted-tree.qza
# qiime tools export \
# --input-path 4_rooted-tree.qza \
# --output-path exported-tree/4_rooted-tree
# 
# #3_unrooted-tree.qza
# qiime tools export \
# --input-path 3_unrooted-tree.qza \
# --output-path exported-tree/3_unrooted-tree


########## betaNTI
library(picante)
library(ape)
library(NST)

setwd("~/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/20250207/20250207_phylo/")


### load otutab and phylotree
# Formatting:
# 1)otutab: row~sample
# 2) phylogree: rooted and Tip_label without $’$.
otu<-  otu_table(ps.filt) %>%t%>% data.frame()
group <- sample_data(ps.filt) %>% data.frame()
group<- group %>% 
  mutate(treatment =factor( paste(Geno,Compartment, sep="_"))) %>%
  select(treatment)
phylo<-read.tree("/Users/hupeiyao/Documents/Research/LabMember/xjw/Publish/MP/20250207/filt/20250207/20250207_phylo/exported-tree/4_rooted-tree/tree.nwk")
phylo$tip.label <- gsub("['‘’]", "", phylo$tip.label)

### pNST
set.seed(123)
pnst <- pNST(comm = otu, tree = phylo, group = group
             , phylo.shuffle = TRUE, taxo.null.model = NULL, 
             pd.wd = tempdir(), abundance.weighted = TRUE, rand = 999, nworker = 4, SES = TRUE, RC = T)

betaMNTD <- pnst$index.pair
head(betaMNTD)
write.csv(betaMNTD, 'betaMNTD.csv', quote = FALSE, row.names= FALSE)

Kitaake_Root <- rownames(subset(group, treatment=='Kitaake_Root'))
betaMNTD_Kitaake_Root <- subset(betaMNTD, name1 %in% Kitaake_Root & name2 %in% Kitaake_Root)
MH63_Root <- rownames(subset(group, treatment=='MH63_Root'))
betaMNTD_MH63_Root <- subset(betaMNTD, name1 %in% MH63_Root & name2 %in% MH63_Root)
Kitaake_Rhizo <- rownames(subset(group, treatment=='Kitaake_Rhizo'))
betaMNTD_Kitaake_Rhizo <- subset(betaMNTD, name1 %in% Kitaake_Rhizo & name2 %in% Kitaake_Rhizo)
MH63_Rhizo <- rownames(subset(group, treatment=='MH63_Rhizo'))
betaMNTD_MH63_Rhizo <- subset(betaMNTD, name1 %in% MH63_Rhizo & name2 %in% MH63_Rhizo)

#|βNTI| < 2, number of occurrences / total number = contribution rate of stochastic factors
nrow(betaMNTD_Kitaake_Root[which(abs(betaMNTD_Kitaake_Root$bNTI.wt)<2), ])/nrow(betaMNTD_Kitaake_Root)  #Control 组
nrow(betaMNTD_MH63_Root[which(abs(betaMNTD_MH63_Root$bNTI.wt)<2), ])/nrow(betaMNTD_MH63_Root)  #Treat 组

#|βNTI| > 2, number of occurrences / total number = contribution rate of deterministic factors
nrow(betaMNTD_Kitaake_Root[which(abs(betaMNTD_Kitaake_Root$bNTI.wt)>2), ])/nrow(betaMNTD_Kitaake_Root)  #Control 组
nrow(betaMNTD_MH63_Root[which(abs(betaMNTD_MH63_Root$bNTI.wt)>2), ])/nrow(betaMNTD_MH63_Root)  #Treat 组

###
library(ggpubr)

betaMNTD_Kitaake_Root$group <- 'Kitaake_Root'
betaMNTD_MH63_Root$group <- 'MH63_Root'
betaMNTD_Kitaake_Rhizo$group <- 'Kitaake_Rhizo'
betaMNTD_MH63_Rhizo$group <- 'MH63_Rhizo'

colorRGB<- c("#795548","#1A237E","#006064",
             "#D4938A","#3F51B5","#00BCD4")

betaMNTD_group <- rbind(betaMNTD_Kitaake_Root, betaMNTD_MH63_Root,betaMNTD_Kitaake_Rhizo, betaMNTD_MH63_Rhizo)
p1<- betaMNTD_group %>% mutate(group = factor(group, levels = c("MH63_Rhizo","MH63_Root","Kitaake_Rhizo","Kitaake_Root"))) %>%
  ggboxplot(x = 'group', y = 'bNTI.wt', color = 'group', fill="group", alpha=.6,width=.5) +
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = list(c('Kitaake_Root', 'MH63_Root'),c('Kitaake_Rhizo', 'MH63_Rhizo'),
                                        c('Kitaake_Root', 'Kitaake_Rhizo'),c('MH63_Root', 'MH63_Rhizo'))) +
  labs(y = 'betaNTI', x="")+geom_jitter(width=.12, aes(color=group))+scale_color_manual(values = colorRGB[c(2,3,5,6)])+scale_fill_manual(values = colorRGB[c(2,3,5,6)])
# ggsave("Figure1.betaNTI.pdf", width=5, height = 5)

p2<- betaMNTD_group %>% mutate(group = factor(group, levels = c("MH63_Rhizo","MH63_Root","Kitaake_Rhizo","Kitaake_Root"))) %>%
  ggboxplot(x = 'group', y = 'RC.bMNTD.wt', color = 'group', fill="group", alpha=.6,width=.5) +
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = list(c('Kitaake_Root', 'MH63_Root'),c('Kitaake_Rhizo', 'MH63_Rhizo'),
                                        c('Kitaake_Root', 'Kitaake_Rhizo'),c('MH63_Root', 'MH63_Rhizo'))) +
  labs(y = 'RCbray', x="")+geom_jitter(width=.12, aes(color=group))+scale_color_manual(values = colorRGB[c(2,3,5,6)])+scale_fill_manual(values = colorRGB[c(2,3,5,6)])


#### composition
library(dplyr)
betaMNTD_group <- betaMNTD_group %>%
  mutate(
    process = case_when(
      bNTI.wt > 2 ~ "Variable selection",
      bNTI.wt < -2 ~ "Homogeneous selection",
      abs(bNTI.wt) < 2 & abs(RC.bMNTD.wt) > 0.95 & RC.bMNTD.wt > 0.95 ~ "Dispersal limitation",
      abs(bNTI.wt) < 2 & abs(RC.bMNTD.wt) > 0.95 & RC.bMNTD.wt < -0.95 ~ "Homogenizing dispersal",
      abs(bNTI.wt) < 2 & abs(RC.bMNTD.wt) < 0.95 ~ "Drift alone",
      TRUE ~ "Unknown"
    )
  )

betaMNTD_group$process <- factor(betaMNTD_group$process, levels = c("Homogeneous selection", "Dispersal limitation", "Homogenizing dispersal", "Drift alone"))

p3<- betaMNTD_group %>%
  mutate(group = factor(group, levels = c("MH63_Rhizo","MH63_Root","Kitaake_Rhizo","Kitaake_Root"))) %>%
  select(group, process) %>%
  count(process, group) %>%
  group_by(group) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  select(group, process, percentage) %>%
  ggplot(aes(x = group, y = percentage, fill = process)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("skyblue1", "deeppink1", "burlywood1", "springgreen3", "grey60")) +
  labs(title = "") +
  xlab("") + ylab("ASV_betaNTI_RCbray\nAssembly process (%)") +
  theme_pubr()


p<- gridExtra::grid.arrange(p1, p2, p3,ncol = 3)
ggsave("Figure1.AssemblyEvaluation.pdf",p, width=15, height=5)

save.image("20250208.deteministic_stomastic_process.RData")

