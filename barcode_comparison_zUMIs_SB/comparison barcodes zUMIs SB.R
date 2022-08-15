### code for template analysis first scRNAseq experiment with SPECTRAL

### files used for example --> mouse unmod human unmod ORF

### to access files --> move files to /home/gmoro/RStudio/SPECTRAL


########
#STANDARD WORKFLOW TO GENERATE UMAPS
########

### clearing environment 

rm(list=ls())

### loading packages 

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scDblFinder)
library(BiocManager)
library(tidyr)
library(gridExtra)
library(stringr)
library(seqinr)
library(Biostrings)
library(argparse)

# table to change barcodes 

cell_labels_2Comm <- readDNAStringSet("~/RStudio/cell_labels_2Comm.fasta")

dna.sequences <- data.frame(ID=names(cell_labels_2Comm),sequences=as.character(cell_labels_2Comm))


### data zUMIs

z_100_100_TSO<- readRDS("~/SPECTRAL/sixth_scRNAseq/100_100_TSO.dgecounts.rds")
z_100_100_dT<- readRDS("~/SPECTRAL/sixth_scRNAseq/100_100_dT.dgecounts.rds")
z_100_100_combined<- readRDS("~/SPECTRAL/sixth_scRNAseq/100_100_combined.dgecounts.rds")

z_100_10_TSO<- readRDS("~/SPECTRAL/sixth_scRNAseq/100_10_TSO.dgecounts.rds")
z_100_10_dT<- readRDS("~/SPECTRAL/sixth_scRNAseq/100_10_dT.dgecounts.rds")
z_100_10_combined<- readRDS("~/SPECTRAL/sixth_scRNAseq/100_10_combined.dgecounts.rds")

z_10_100_TSO<- readRDS("~/SPECTRAL/sixth_scRNAseq/10_100_TSO.dgecounts.rds")
z_10_100_dT<- readRDS("~/SPECTRAL/sixth_scRNAseq/10_100_dT.dgecounts.rds")
z_10_100_combined<- readRDS("~/SPECTRAL/sixth_scRNAseq/10_100_combined.dgecounts.rds")

z_unmod_TSO<- readRDS("~/SPECTRAL/sixth_scRNAseq/unmod_TSO.dgecounts.rds")
z_unmod_dT<- readRDS("~/SPECTRAL/sixth_scRNAseq/unmod_dT.dgecounts.rds")
z_unmod_combined<- readRDS("~/SPECTRAL/sixth_scRNAseq/unmod_combined.dgecounts.rds")

### subsetting data zUMIs

z_100_100_TSO<-z_100_100_TSO$umicount$exon$all
length(colnames(z_100_100_TSO))

z_100_100_dT<-z_100_100_dT$umicount$exon$all
length(colnames(z_100_100_dT))

z_100_100_combined<-z_100_100_combined$umicount$exon$all
length(colnames(z_100_100_combined))



z_100_10_TSO<-z_100_10_TSO$umicount$exon$all
length(colnames(z_100_10_TSO))

z_100_10_dT<-z_100_10_dT$umicount$exon$all
length(colnames(z_100_10_dT))

z_100_10_combined<-z_100_10_combined$umicount$exon$all
length(colnames(z_100_10_combined))


z_10_100_TSO<-z_10_100_TSO$umicount$exon$all
length(colnames(z_10_100_TSO))

z_10_100_dT<-z_10_100_dT$umicount$exon$all
length(colnames(z_10_100_dT))

z_10_100_combined<-z_10_100_combined$umicount$exon$all
length(colnames(z_10_100_combined))


z_unmod_TSO<-z_unmod_TSO$umicount$exon$all
length(colnames(z_unmod_TSO))

z_unmod_dT<-z_unmod_dT$umicount$exon$all
length(colnames(z_unmod_dT))

z_unmod_combined<-z_unmod_combined$umicount$exon$all
length(colnames(z_unmod_combined))


# changing gene names

gene_names_z_100_100_TSO <- read.delim("~/SPECTRAL/sixth_scRNAseq/100_100_TSO.gene_names.txt")

rownames(z_100_100_TSO)<-
  dplyr::recode(
    rownames(z_100_100_TSO),
    !!!setNames(as.character(gene_names_z_100_100_TSO$gene_name),gene_names_z_100_100_TSO$gene_id)
  )

gene_names_z_100_100_dT <- read.delim("~/SPECTRAL/sixth_scRNAseq/100_100_dT.gene_names.txt")

rownames(z_100_100_dT)<-
  dplyr::recode(
    rownames(z_100_100_dT ),
    !!!setNames(as.character(gene_names_z_100_100_dT $gene_name),gene_names_z_100_100_dT$gene_id)
  )

gene_names_z_100_100_combined <- read.delim("~/SPECTRAL/sixth_scRNAseq/100_100_combined.gene_names.txt")

rownames(z_100_100_combined)<-
  dplyr::recode(
    rownames(z_100_100_combined),
    !!!setNames(as.character(gene_names_z_100_100_combined$gene_name),gene_names_z_100_100_combined$gene_id)
  )


gene_names_z_100_10_TSO <- read.delim("~/SPECTRAL/sixth_scRNAseq/100_10_TSO.gene_names.txt")

rownames(z_100_10_TSO)<-
  dplyr::recode(
    rownames(z_100_10_TSO),
    !!!setNames(as.character(gene_names_z_100_10_TSO$gene_name),gene_names_z_100_10_TSO$gene_id)
  )

gene_names_z_100_10_dT <- read.delim("~/SPECTRAL/sixth_scRNAseq/100_10_dT.gene_names.txt")

rownames(z_100_10_dT)<-
  dplyr::recode(
    rownames(z_100_10_dT ),
    !!!setNames(as.character(gene_names_z_100_10_dT $gene_name),gene_names_z_100_10_dT$gene_id)
  )

gene_names_z_100_10_combined <- read.delim("~/SPECTRAL/sixth_scRNAseq/100_10_combined.gene_names.txt")

rownames(z_100_10_combined)<-
  dplyr::recode(
    rownames(z_100_10_combined),
    !!!setNames(as.character(gene_names_z_100_10_combined$gene_name),gene_names_z_100_10_combined$gene_id)
  )



gene_names_z_10_100_TSO <- read.delim("~/SPECTRAL/sixth_scRNAseq/10_100_TSO.gene_names.txt")

rownames(z_10_100_TSO)<-
  dplyr::recode(
    rownames(z_10_100_TSO),
    !!!setNames(as.character(gene_names_z_10_100_TSO$gene_name),gene_names_z_10_100_TSO$gene_id)
  )

gene_names_z_10_100_dT <- read.delim("~/SPECTRAL/sixth_scRNAseq/10_100_dT.gene_names.txt")

rownames(z_10_100_dT)<-
  dplyr::recode(
    rownames(z_10_100_dT ),
    !!!setNames(as.character(gene_names_z_10_100_dT $gene_name),gene_names_z_10_100_dT$gene_id)
  )

gene_names_z_10_100_combined <- read.delim("~/SPECTRAL/sixth_scRNAseq/10_100_combined.gene_names.txt")

rownames(z_10_100_combined)<-
  dplyr::recode(
    rownames(z_10_100_combined),
    !!!setNames(as.character(gene_names_z_10_100_combined$gene_name),gene_names_z_10_100_combined$gene_id)
  )



gene_names_z_unmod_TSO <- read.delim("~/SPECTRAL/sixth_scRNAseq/unmod_TSO.gene_names.txt")

rownames(z_unmod_TSO)<-
  dplyr::recode(
    rownames(z_unmod_TSO),
    !!!setNames(as.character(gene_names_z_unmod_TSO$gene_name),gene_names_z_unmod_TSO$gene_id)
  )

gene_names_z_unmod_dT <- read.delim("~/SPECTRAL/sixth_scRNAseq/unmod_dT.gene_names.txt")

rownames(z_unmod_dT)<-
  dplyr::recode(
    rownames(z_unmod_dT ),
    !!!setNames(as.character(gene_names_z_unmod_dT $gene_name),gene_names_z_unmod_dT$gene_id)
  )

gene_names_z_unmod_combined <- read.delim("~/SPECTRAL/sixth_scRNAseq/unmod_combined.gene_names.txt")

rownames(z_unmod_combined)<-
  dplyr::recode(
    rownames(z_unmod_combined),
    !!!setNames(as.character(gene_names_z_unmod_combined$gene_name),gene_names_z_unmod_combined$gene_id)
  )



### data SB

sb_100_100_dT<-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/grouped_GTGA_GACA_Tn_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_100_100_TSO<-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/grouped_AATG_CCAC_TSO_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_100_100_combined<-read.csv("~/RStudio/SB/100_mod_100_primer_no_cap_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

sb_unmod_TSO <-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/unmod_TSO_no_cap_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_unmod_dT <-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/unmod_no_cap_Tn_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_unmod_combined<-read.csv("~/RStudio/SB/unmod_no_cap_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

sb_10_100_TSO  <-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/no_cap_10%_100uM_TSO_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_10_100_dT <-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/no_cap_10%_100uM_Tn_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_10_100_combined<-read.csv("~/RStudio/SB/10_mod_100_primer_no_cap_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

sb_100_10_TSO <-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/no_cap_100%_mod_10uM_TSO_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_100_10_dT <-read.csv("~/RStudio/SB/100_mod_100_primer_sep_reads/no_cap_100%_mod_10uM_Tn_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")
sb_100_10_combined<-read.csv("~/RStudio/SB/final_100_mod_10_primer_RSEC_MolsPerCell.csv", row.names=1, comment.char="#")

##### converting cell barcodes sb to zUMIs

cb1 <- read.table("/home/gmoro/RStudio/whitelist/CLS1")$V1
cb2 <- read.table("/home/gmoro/RStudio/whitelist/CLS2")$V1
cb3 <- read.table("/home/gmoro/RStudio/whitelist/CLS3")$V1

### refined code

tests =c(130633,43851,61929,138720,139200)

translate_id <- function(x, cb1, cb2, cb3) {
  cl1_int <- floor(x / (96**2)) + 1 + 1
  cl2_int <- (floor(x %% (96**2)) / 96) + 1 + 1
  cl3_int <- (floor(x %% (96**2)) %% 96) + 1
  if (cl3_int == 1){
    cl3_int <- 97 + 1
  }
  cl1 <- cb1[cl1_int -1] 
  cl2 <- cb2[cl2_int -1]
  cl3 <- cb3[cl3_int -1]
  cl <- paste(c(cl1, cl2, cl3), collapse = '-')
}
  

for (idx in tests){
  print(translate_id(x = idx, cb1, cb2 = cb2, cb3 = cb3))
}

### test with actual barcodes

rownames(sb_100_100_combined)

barcodes_z<-c()

for (i in rownames(sb_100_100_combined)) {
  i <- as.integer(i)
  print(i)
  cl1_int <- floor(i / (96**2)) + 1 + 1
  cl2_int <- (floor(i %% (96**2)) / 96) + 1 + 1
  cl3_int <- (floor(i %% (96**2)) %% 96) + 1
  if (cl3_int == 1){
    cl3_int <- 97 + 1
  }
  cl1 <- cb1[cl1_int -1] 
  cl2 <- cb2[cl2_int -1]
  cl3 <- cb3[cl3_int -1]
  cl <- paste(c(cl1, cl2, cl3), collapse = '')
  print(cl)
  barcodes_z <- c(barcodes_z,cl)
}

rownames(sb_100_100_combined) <- barcodes_z

sb_100_100_combined

length(intersect(rownames(sb_100_100_combined),colnames(z_100_100_combined))) # 16676

length(rownames(sb_100_100_combined)) # 19844

length(colnames(z_100_100_combined)) # 33880



