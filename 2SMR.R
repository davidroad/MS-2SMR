####2 sample MR for eQTL and selected complex trait.
###Date: 3/24/2021
###Package
###Install and load package
#BiocManager::install("devtools")
#install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("MRCIEU/MRInstruments", force = T)

###loading package
library("GOstats")
library("RDAVIDWebService")
library("org.Hs.eg.db")
library("clusterProfiler")
library(dplyr)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(liftOver)
library(stringr)
library(plyr)
library(readxl)

###Load and organize data
getwd()
anno <- read.delim("GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")

####input selected tissues.
Tissue <-  c("Brain Frontal Cortex BA9","Brain Hippocampus","Brain Hypothalamus",
             "Brain Nucleus accumbens basal ganglia","Brain Putamen basal ganglia",
             "Brain Spinal cord cervical C1","Brain Substantia nigra","Colon Sigmoid",
             "Heart Left Ventricle","Stomach","Spleen", 
             "Whole Blood","Brain Cerebellum","Artery Coronary","Brain Amygdala",
             "Brain Anterior cingulate cortex BA24",
             "Brain Caudate basal ganglia","Brain Cerebellar Hemisphere","Brain Cortex")


####Conducting 2 sample MR for selected tissues eqtl
for (i in Tissue){
  cat("Starting 2SMR for", i,"\n")
  gc()
  addr <- paste("./Tissues_GTEx/",
                i,".txt", 
                sep = "")
  
  ###Read data and add header
  exp_data <- read.delim(addr,header = F)
  colnames(exp_data) <- c("ENSEMBL","var_id","dist_tss","ma_samples","ma_count",
                          "eaf","pval","beta","se")
  
  cat(length(unique(exp_data$ENSEMBL)),"genes have probes passed p<0.0001","\n")
  
  ###annotate variant id to snp id
  exp_data <- merge.data.frame(exp_data,anno, by.x = "var_id",by.y = "variant_id")
  exp_data$snp_position <- paste(str_split_fixed(exp_data$chr,"r",n =2)[,2],":",exp_data$variant_pos,sep = "")
  
  exp_data$ENSEMBL <- str_split_fixed(exp_data$ENSEMBL, "\\.",n = 2)[,1]
  
  tb <- unique(exp_data$ENSEMBL)
  tb = bitr(tb, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
  
  exp_data <- merge.data.frame(exp_data, tb, by.x = "ENSEMBL", by.y = "ENSEMBL")
  
  
  final <- exp_data[,c("SYMBOL","rs_id_dbSNP151_GRCh38p7",
                       "snp_position","alt","ref","beta","se","pval","eaf")]
  
  colnames(final) <- c("gene_name","SNP","snp_position",
                       "effect_allele","other_allele","beta","se","pval","eaf")
  
  ###Tissue
  final$tissue <- i
  final <- final[complete.cases(final),]
  
  ###Exposure
  ###Formating based on 2SMR pacakege
  cat("formating exposure:", i,"...")
  eqtl_exp_dat <- format_gtex_eqtl(subset(final))
  eqtl_exp_dat <- eqtl_exp_dat[complete.cases(eqtl_exp_dat),]
  
  ###clumping
  cat("clumping...may take longer","\n")
  eqtl_exp_dat <- clump_data(eqtl_exp_dat)
  
  ###Outcome
  cat("finding outcome...","\n")
  out_dat <- extract_outcome_data(
    snps = eqtl_exp_dat$SNP,
    outcomes = 'ieu-b-18'
  )
  
  ###Harmoise Data
  cat("harmonising...","\n")
  dat <- harmonise_data(
    exposure_dat = eqtl_exp_dat,
    outcome_dat = out_dat
  )
  
  ###Run MR analysis
  cat("runing 2SMR...","\n")
  res <- mr(dat)
   
  ###OUTPUT 2SMR result
  cat("writing table...","\n")
  
  out1 <- paste("./2SMR results/",
                       i,"_2SMR.csv",
                       sep = "")
  
  write.csv(res,out1)
  
  ###Finding significant result
  res$fdr <- p.adjust(res$pval, method = "BH", n = length(res$pval))
  res$gene <- str_split_fixed(res$exposure, " ", n = 2)[,1]
  res_sig <- res[res$fdr < 0.05,]
  res_sig <- res_sig[,c("gene","outcome","method","nsnp", "b","se", "pval","fdr")]
  
  ###OUTPUT 2SMR SIG result
  out2 <- paste("./Sig/",
                i,"_2SMR_sig.csv",
                sep = "")
  write.csv(res_sig, out2)
  
  
  
  #############for manhattan plot
  eqtl_exp_dat$gene <- str_split_fixed(eqtl_exp_dat$exposure, " ", n = 2)[,1]
  
  pre <- merge(final,eqtl_exp_dat,
               by.x = c("SNP","gene_name"),
               by.y = c("SNP","gene"))
  
  cob <- merge(res, pre, by = "id.exposure")
  cob$chr <- as.numeric(str_split_fixed(cob$snp_position,":",n = 2)[,1])
  cob$bp <- as.numeric(str_split_fixed(cob$snp_position,":",n = 2)[,2])
  
  ###Loading Bayesian risk genes
  risk_gene <- read.delim("MS_Risk_GeneList_methy.txt")
  rlist <- risk_gene$x
  
  ###
  inc <- match(cob$gene_name,rlist)
  cob$highlight <- 0
  cob[cob$gene_name %in% rlist & cob$fdr < 0.05,]$highlight <- 1 
  
  
  addr <- paste("./Manhattan/",
                i,"_2SMR_Manhattan.csv",
                sep = "")
  
  write.csv(cob,file = addr)
  
  
  
  cat(i,"Done")
}



for (i in Tissue){
  addr <- paste("./2SMR results/",i,"_2SMR.txt",sep = "")
  df <- read.delim(addr, sep = ",", row.names = 1)
  df$FDR <- p.adjust(df$pval, method = "BH",n = length(df$pval))
  
  addr1 <- paste("./2SMR results/",i,"_2SMR.csv",sep = "")
  write.csv(df, addr1)
}


