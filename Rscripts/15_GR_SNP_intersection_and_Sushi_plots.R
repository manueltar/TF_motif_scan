
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

select_motifs = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform SNP ----
  
  SNP = opt$SNP
  
  cat("SNP_\n")
  cat(sprintf(as.character(SNP)))
  cat("\n")
  
  chr<-paste('chr',gsub(":.+$","",SNP),sep='')
  START<-as.integer(gsub("^[^:]+:","",SNP))
  
  
  cat("SNP\n")
  cat(sprintf(as.character(chr)))
  cat("\t")
  cat(sprintf(as.character(START)))
  cat("\n")
  
  gr_SNP <- GRanges(
    seqnames = as.character(chr),
    name2='SNP',
    name3='SNP',
    ranges=IRanges(
      start=START,
      names='SNP',
      end=START,
      strand="*"))
  
  cat("gr_SNP_0\n")
  cat(str(gr_SNP))
  cat("\n")
  
  #### READ and transform ROI ----
  
  ROI = data.frame(unlist(strsplit(opt$ROI, split=",")), stringsAsFactors = F)
  
  colnames(ROI)<-'interval'
  
  cat("ROI_0\n")
  cat(str(ROI))
  cat("\n")
  
  chr<-paste('chr',gsub(":.+$","",ROI$interval),sep='')
  START<-gsub("-.+$","",ROI$interval)
  START<-as.integer(gsub("^[^:]+:","",START))
  END<-gsub("^[^-]+-","",ROI$interval)
  END<-as.integer(gsub("__.+$","",END))
  names_ROI<-gsub("^[^__]+__","",ROI$interval)
  
  START_MASTER<-START[1]
  chr_MASTER<-unique(chr)
  
  
  cat("ROI\n")
  cat(sprintf(as.character(chr)))
  cat("\t")
  cat(sprintf(as.character(START)))
  cat("\t")
  cat(sprintf(as.character(END)))
  cat("\n")
  cat(sprintf(as.character(names_ROI)))
  cat("\n")
  
  gr_ROI <- GRanges(
    seqnames = as.character(chr),
    name2='ROI',
    name3='ROI',
    ranges=IRanges(
      start=START,
      names=names_ROI,
      end=END,
      strand="*"))
  
  # cat("gr_ROI_0\n")
  # cat(str(gr_ROI))
  # cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ TF_collapsed_EXP_crossed ----
  
  TF_collapsed_EXP_crossed<-as.data.frame(fread(file=opt$TF_collapsed_EXP_crossed, sep="\t", header=T), stringsAsFactors=F)
  
  # Discard motifs only found in CIS_BP
  
  TF_collapsed_EXP_crossed<-unique(TF_collapsed_EXP_crossed[which(TF_collapsed_EXP_crossed$string_source != 'CIS_BP'),])
  
  cat("TF_collapsed_EXP_crossed_0\n")
  cat(str(TF_collapsed_EXP_crossed))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_crossed$Symbol)))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_crossed$ensembl_gene_id)))
  cat("\n")
  
  TF_collapsed_EXP_crossed.dt<-data.table(TF_collapsed_EXP_crossed, key=c("Symbol","ensembl_gene_id"))
  
  TF_collapsed_EXP_crossed_expanded<-as.data.frame(TF_collapsed_EXP_crossed.dt[,.(source=unlist(strsplit(string_source, split=';'))), by=key(TF_collapsed_EXP_crossed.dt)], stringsAsFactors=F)
  
  
  cat("TF_collapsed_EXP_crossed_expanded_0\n")
  cat(str(TF_collapsed_EXP_crossed_expanded))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_crossed_expanded$Symbol)))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_crossed_expanded$ensembl_gene_id)))
  cat("\n")
  
  #### READ TF_collapsed_EXP_and_DE_crossed ----
  
  TF_collapsed_EXP_and_DE_crossed<-as.data.frame(fread(file=opt$TF_collapsed_EXP_and_DE_crossed, sep="\t", header=T), stringsAsFactors=F)
  
  # Discard motifs only found in CIS_BP
  
  TF_collapsed_EXP_and_DE_crossed<-unique(TF_collapsed_EXP_and_DE_crossed[which(TF_collapsed_EXP_and_DE_crossed$string_source != 'CIS_BP'),])
  
  cat("TF_collapsed_EXP_and_DE_crossed_0\n")
  cat(str(TF_collapsed_EXP_and_DE_crossed))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_and_DE_crossed$Symbol)))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_and_DE_crossed$ensembl_gene_id)))
  cat("\n")
  
  TF_collapsed_EXP_and_DE_crossed.dt<-data.table(TF_collapsed_EXP_and_DE_crossed, key=c("Symbol","ensembl_gene_id"))
  
  TF_collapsed_EXP_and_DE_crossed_expanded<-as.data.frame(TF_collapsed_EXP_and_DE_crossed.dt[,.(source=unlist(strsplit(string_source, split=';'))), by=key(TF_collapsed_EXP_and_DE_crossed.dt)], stringsAsFactors=F)
  
  cat("TF_collapsed_EXP_and_DE_crossed_expanded_0\n")
  cat(str(TF_collapsed_EXP_and_DE_crossed_expanded))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_and_DE_crossed_expanded$Symbol)))
  cat("\n")
  cat(str(unique(TF_collapsed_EXP_and_DE_crossed_expanded$ensembl_gene_id)))
  cat("\n")
  
 
  #### READ input_file_HOMER ----
  
  input_file_HOMER<-as.data.frame(fread(file=opt$input_file_HOMER, sep="\t", header=F, skip=5), stringsAsFactors=F)
  colnames(input_file_HOMER)<-c('query_region','start','end','Motif_ID','score','strand')
  
  cat("input_file_HOMER_0\n")
  cat(str(input_file_HOMER))
  cat("\n")
  cat(str(unique(input_file_HOMER$Motif_ID)))
  cat("\n")
  
  input_file_HOMER$Symbol<-input_file_HOMER$Motif_ID
  
  input_file_HOMER$Symbol<-gsub("/.+$","",input_file_HOMER$Symbol)
  input_file_HOMER$Symbol<-gsub("^[^_]+_","",input_file_HOMER$Symbol)
  input_file_HOMER$Symbol<-gsub("\\(.+$","",input_file_HOMER$Symbol)
  input_file_HOMER$Symbol<-gsub("-halfsite","",input_file_HOMER$Symbol)
  
  input_file_HOMER$source<-'HOMER'
  
  
  cat("input_file_HOMER_3\n")
  cat(str(input_file_HOMER))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_HOMER$Symbol))))
  cat("\n")
  
  
  input_file_HOMER_sel_TF_EXP<-merge(input_file_HOMER,
                                TF_collapsed_EXP_crossed_expanded,
                                by=c("Symbol","source"))
  
  
  cat("input_file_HOMER_sel_TF_EXP_0\n")
  cat(str(input_file_HOMER_sel_TF_EXP))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_HOMER_sel_TF_EXP$Symbol))))
  cat("\n")
  
  input_file_HOMER_sel_TF_EXP_DE<-merge(input_file_HOMER,
                                     TF_collapsed_EXP_and_DE_crossed_expanded,
                                     by=c("Symbol","source"))
  
  
  cat("input_file_HOMER_sel_TF_EXP_DE_0\n")
  cat(str(input_file_HOMER_sel_TF_EXP_DE))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_HOMER_sel_TF_EXP_DE$Symbol))))
  cat("\n")
  
  
  #### READ input_file_JASPAR ----
  
  input_file_JASPAR<-as.data.frame(fread(file=opt$input_file_JASPAR, sep="\t", header=F, skip=5), stringsAsFactors=F)
  colnames(input_file_JASPAR)<-c('query_region','start','end','Motif_ID','score','strand')
  
  cat("input_file_JASPAR_0\n")
  cat(str(input_file_JASPAR))
  cat("\n")
  cat(str(unique(input_file_JASPAR$Motif_ID)))
  cat("\n")
  
  input_file_JASPAR$Symbol<-input_file_JASPAR$Motif_ID
  
  input_file_JASPAR$Symbol<-gsub("^[^_]+_","",input_file_JASPAR$Symbol)
  input_file_JASPAR$Symbol<-gsub("\\(.+$","",input_file_JASPAR$Symbol)
  
    
  # cat("input_file_JASPAR_2\n")
  # cat(str(input_file_JASPAR))
  # cat("\n")
  # cat(sprintf(as.character(unique(input_file_JASPAR$Symbol))))
  # cat("\n")
  
  input_file_JASPAR.dt<-data.table(input_file_JASPAR, key=c('query_region','start','end','Motif_ID','score','strand'))
  
  input_file_JASPAR_DEF<-as.data.frame(input_file_JASPAR.dt[,.(Symbol=unlist(strsplit(Symbol, split='::'))), by = key(input_file_JASPAR.dt)], stringsAsFactors=F)
  input_file_JASPAR_DEF$source<-'JASPAR'
  
  
  cat("input_file_JASPAR_DEF_2\n")
  cat(str(input_file_JASPAR_DEF))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_JASPAR_DEF$Symbol))))
  cat("\n")
  
  input_file_JASPAR_DEF_sel_TF_EXP<-merge(input_file_JASPAR_DEF,
                                     TF_collapsed_EXP_crossed_expanded,
                                     by=c("Symbol","source"))
  
  
  cat("input_file_JASPAR_DEF_sel_TF_EXP_0\n")
  cat(str(input_file_JASPAR_DEF_sel_TF_EXP))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_JASPAR_DEF_sel_TF_EXP$Symbol))))
  cat("\n")
  
  input_file_JASPAR_DEF_sel_TF_EXP_DE<-merge(input_file_JASPAR_DEF,
                                        TF_collapsed_EXP_and_DE_crossed_expanded,
                                        by=c("Symbol","source"))
  
  
  cat("input_file_JASPAR_DEF_sel_TF_EXP_DE_0\n")
  cat(str(input_file_JASPAR_DEF_sel_TF_EXP_DE))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_JASPAR_DEF_sel_TF_EXP_DE$Symbol))))
  cat("\n")
  
  #### READ input_file_CIS_BP ----
  
  input_file_CIS_BP<-as.data.frame(fread(file=opt$input_file_CIS_BP, sep="\t", header=T), stringsAsFactors=F)

  input_file_CIS_BP$Symbol<-input_file_CIS_BP$TF_Name
  
  input_file_CIS_BP$source<-'CIS_BP'
  
  
  cat("input_file_CIS_BP_0\n")
  cat(str(input_file_CIS_BP))
  cat("\n")
  cat(str(unique(input_file_CIS_BP$Motif_ID)))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_CIS_BP$Symbol))))
  cat("\n")
  
  input_file_CIS_BP_sel_TF_EXP<-merge(input_file_CIS_BP,
                                     TF_collapsed_EXP_crossed_expanded,
                                     by=c("Symbol","source"))
  
  
  cat("input_file_CIS_BP_sel_TF_EXP_0\n")
  cat(str(input_file_CIS_BP_sel_TF_EXP))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_CIS_BP_sel_TF_EXP$Symbol))))
  cat("\n")
  
  input_file_CIS_BP_sel_TF_EXP_DE<-merge(input_file_CIS_BP,
                                        TF_collapsed_EXP_and_DE_crossed_expanded,
                                        by=c("Symbol","source"))
  
  
  cat("input_file_CIS_BP_sel_TF_EXP_DE_0\n")
  cat(str(input_file_CIS_BP_sel_TF_EXP_DE))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_CIS_BP_sel_TF_EXP_DE$Symbol))))
  cat("\n")
  
  #### Get the TF expanded motifs together ----
  
  
  
  gr_HOMER_TF_EXP <- GRanges(
    seqnames = as.character(chr_MASTER),
    name2='HOMER',
    name3=input_file_HOMER_sel_TF_EXP$Motif_ID,
    name4=input_file_HOMER_sel_TF_EXP$score,
    name5=input_file_HOMER_sel_TF_EXP$query_region,
    strand=input_file_HOMER_sel_TF_EXP$strand,
    ranges=IRanges(
      start=START_MASTER+input_file_HOMER_sel_TF_EXP$start,
      end=START_MASTER+input_file_HOMER_sel_TF_EXP$end,
      names=input_file_HOMER_sel_TF_EXP$Symbol))
  
  # cat("gr_HOMER_TF_EXP_0\n")
  # cat(str(gr_HOMER_TF_EXP))
  # cat("\n")
  
  
  
  gr_JASPAR_TF_EXP <- GRanges(
    seqnames = as.character(chr_MASTER),
    strand=input_file_JASPAR_DEF_sel_TF_EXP$strand,
    name2='JASPAR',
    name3=input_file_JASPAR_DEF_sel_TF_EXP$Motif_ID,
    name4=input_file_JASPAR_DEF_sel_TF_EXP$score,
    name5=input_file_JASPAR_DEF_sel_TF_EXP$query_region,
    ranges=IRanges(
      start=START_MASTER+input_file_JASPAR_DEF_sel_TF_EXP$start,
      end=START_MASTER+input_file_JASPAR_DEF_sel_TF_EXP$end,
      names=input_file_JASPAR_DEF_sel_TF_EXP$Symbol))
  
  # cat("gr_JASPAR_TF_EXP_0\n")
  # cat(str(gr_JASPAR_TF_EXP))
  # cat("\n")
  
  gr_CIS_BP_TF_EXP <- unique(GRanges(
    seqnames = as.character(chr_MASTER),
    strand=input_file_CIS_BP_sel_TF_EXP$strand,
    name2='CIS_BP',
    name3=input_file_CIS_BP_sel_TF_EXP$Motif_ID,
    name4=input_file_CIS_BP_sel_TF_EXP$score,
    name5=input_file_CIS_BP_sel_TF_EXP$query_region,
    ranges=IRanges(
      start=START_MASTER+input_file_CIS_BP_sel_TF_EXP$start,
      end=START_MASTER+input_file_CIS_BP_sel_TF_EXP$end,
      names=input_file_CIS_BP_sel_TF_EXP$Symbol)))
  
  # cat("gr_CIS_BP_TF_EXP_0\n")
  # cat(str(gr_CIS_BP_TF_EXP))
  # cat("\n")
  
  gr_TF_EXP<-append(gr_HOMER_TF_EXP,gr_JASPAR_TF_EXP)
  
  # cat("gr_TF_EXP_0\n")
  # cat(str(gr_TF_EXP))
  # cat("\n")
  
  gr_TF_EXP<-append(gr_TF_EXP,gr_CIS_BP_TF_EXP)
  
  cat("gr_TF_EXP_1\n")
  cat(str(gr_TF_EXP))
  cat("\n")
  
  #### find overlap with SNP & ROI---
  
  m <- findOverlaps(gr_SNP,gr_TF_EXP)
  
  # cat("m\n")
  # cat(str(m))
  # cat("\n")
  
  subjectHits_SNP<-subjectHits(m)
  
  cat("subjectHits_SNP\n")
  cat(str(subjectHits_SNP))
  cat("\n")
  
  m <- findOverlaps(gr_ROI,gr_TF_EXP)
  
  # cat("m\n")
  # cat(str(m))
  # cat("\n")
  
  subjectHits_ROI<-subjectHits(m)
  
  cat("subjectHits_ROI\n")
  cat(str(subjectHits_ROI))
  cat("\n")
  
  df_SNP <- data.frame(query_region=NA,
                   chr=seqnames(gr_SNP),
                   start=start(gr_SNP),
                   end=end(gr_SNP),
                   strand=c(rep('*', length(gr_SNP))),
                   feature=c(rep('SNP', length(gr_SNP))),
                   source=gr_SNP$name2,
                   Motif_ID=gr_SNP$name3,
                   score=c(rep(0, length(gr_SNP))),
                   Symbol=names(gr_SNP),
                   Intersect_SNP=c(rep('NA', length(gr_SNP))),
                   Intersect_ROI=c(rep('NA', length(gr_SNP))), stringsAsFactors = F)
  # cat("df_SNP_0\n")
  # cat(str(df_SNP))
  # cat("\n")
  
  df_ROI <- data.frame(query_region=NA,
                       chr=seqnames(gr_ROI),
                       start=start(gr_ROI),
                       end=end(gr_ROI),
                       strand=c(rep('*', length(gr_ROI))),
                       feature=names(gr_ROI),
                       source=gr_ROI$name2,
                       Motif_ID=gr_ROI$name3,
                       score=c(rep(0, length(gr_ROI))),
                       Symbol=names(gr_ROI),
                       Intersect_SNP=c(rep('NA', length(gr_ROI))),
                       Intersect_ROI=c(rep('NA', length(gr_ROI))), stringsAsFactors = F)
  
  # cat("df_ROI_0\n")
  # cat(str(df_ROI))
  # cat("\n")
  
  chr_check<-as.character(seqnames(gr_TF_EXP))
  
  cat("chr_check_0\n")
  cat(str(chr_check))
  cat("\n")
  
  STARTS_check<-start(gr_TF_EXP)
  
  cat("STARTS_check_0\n")
  cat(str(STARTS_check))
  cat("\n")
  
  ENDS_check<-end(gr_TF_EXP)
  
  cat("ENDS_check_0\n")
  cat(str(ENDS_check))
  cat("\n")
  
  STRANDS_check<-as.character(strand(gr_TF_EXP))
  
  cat("STRANDS_check_0\n")
  cat(str(STRANDS_check))
  cat("\n")
  
 
  df <- data.frame(query_region=as.character(gr_TF_EXP$name5),
                   chr=as.character(seqnames(gr_TF_EXP)),
                   start=as.integer(start(gr_TF_EXP)),
                   end=as.integer(end(gr_TF_EXP)),
                   strand=as.character(strand(gr_TF_EXP)),
                   feature=c(rep('TF_Motif', length(as.character(seqnames(gr_TF_EXP))))),
                   source=as.character(gr_TF_EXP$name2),
                   Motif_ID=as.character(gr_TF_EXP$name3),
                   score=as.numeric(gr_TF_EXP$name4),
                   Symbol=names(gr_TF_EXP),
                   Intersect_SNP=c(rep('NO', length(as.character(seqnames(gr_TF_EXP))))),
                   Intersect_ROI=c(rep('NO', length(as.character(seqnames(gr_TF_EXP))))), stringsAsFactors = F)
  
  cat("df_0\n")
  cat(str(df))
  cat("\n")
  
  df$Intersect_SNP[subjectHits_SNP]<-'YES'
  df$Intersect_ROI[subjectHits_ROI]<-'YES'
  
  cat("df_1\n")
  cat(str(df))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_ROI))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_ROI)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_SNP))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_SNP)))))
  cat("\n")
  
  check<-df[which(df$Intersect_SNP == 'YES'),]
  
  cat("check\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(unique(check$Symbol))))
  cat("\n")
  
  
  df<-rbind(df_SNP,df_ROI,df)
  
  df$analysis<-"EXP"
  
  cat("df_1\n")
  cat(str(df))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_ROI))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_ROI)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_SNP))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_SNP)))))
  cat("\n")
  
  df_EXP<-df
  
  
  #### Get the TF EXP+DE motifs together ----
  
  gr_HOMER_TF_EXP_DE <- GRanges(
    seqnames = as.character(chr_MASTER),
    strand=input_file_HOMER_sel_TF_EXP_DE$strand,
    name2='HOMER',
    name3=input_file_HOMER_sel_TF_EXP_DE$Motif_ID,
    name4=input_file_HOMER_sel_TF_EXP_DE$score,
    name5=input_file_HOMER_sel_TF_EXP_DE$query_region,
    ranges=IRanges(
      start=START_MASTER+input_file_HOMER_sel_TF_EXP_DE$start,
      end=START_MASTER+input_file_HOMER_sel_TF_EXP_DE$end,
      names=input_file_HOMER_sel_TF_EXP_DE$Symbol))
  
  cat("gr_HOMER_TF_EXP_DE_0\n")
  cat(str(gr_HOMER_TF_EXP_DE))
  cat("\n")
  
  
  gr_JASPAR_TF_EXP_DE <- GRanges(
    seqnames = as.character(chr_MASTER),
    strand=input_file_JASPAR_DEF_sel_TF_EXP_DE$strand,
    name2='JASPAR',
    name3=input_file_JASPAR_DEF_sel_TF_EXP_DE$Motif_ID,
    name4=input_file_JASPAR_DEF_sel_TF_EXP_DE$score,
    name5=input_file_JASPAR_DEF_sel_TF_EXP_DE$query_region,
    ranges=IRanges(
      start=START_MASTER+input_file_JASPAR_DEF_sel_TF_EXP_DE$start,
      end=START_MASTER+input_file_JASPAR_DEF_sel_TF_EXP_DE$end,
      names=input_file_JASPAR_DEF_sel_TF_EXP_DE$Symbol))
  
  cat("gr_JASPAR_TF_EXP_DE_0\n")
  cat(str(gr_JASPAR_TF_EXP_DE))
  cat("\n")
  
  gr_CIS_BP_TF_EXP_DE <- unique(GRanges(
    seqnames = as.character(chr_MASTER),
    strand=input_file_CIS_BP_sel_TF_EXP_DE$strand,
    name2='CIS_BP',
    name3=input_file_CIS_BP_sel_TF_EXP_DE$Motif_ID,
    name4=input_file_CIS_BP_sel_TF_EXP_DE$score,
    name5=input_file_CIS_BP_sel_TF_EXP_DE$query_region,
    ranges=IRanges(
      start=START_MASTER+input_file_CIS_BP_sel_TF_EXP_DE$start,
      end=START_MASTER+input_file_CIS_BP_sel_TF_EXP_DE$end,
      names=input_file_CIS_BP_sel_TF_EXP_DE$Symbol)))
  
  cat("gr_CIS_BP_TF_EXP_DE_0\n")
  cat(str(gr_CIS_BP_TF_EXP_DE))
  cat("\n")
  
  gr_TF_EXP_DE<-append(gr_HOMER_TF_EXP_DE,gr_JASPAR_TF_EXP_DE)
  
  # cat("gr_TF_EXP_DE_0\n")
  # cat(str(gr_TF_EXP_DE))
  # cat("\n")
  
  gr_TF_EXP_DE<-append(gr_TF_EXP_DE,gr_CIS_BP_TF_EXP_DE)
  
  cat("gr_TF_EXP_DE_1\n")
  cat(str(gr_TF_EXP_DE))
  cat("\n")
  
  #### find overlap with SNP & ROI---
  
  m <- findOverlaps(gr_SNP,gr_TF_EXP_DE)
  
  # cat("m\n")
  # cat(str(m))
  # cat("\n")
  
  subjectHits_SNP<-subjectHits(m)
  
  cat("subjectHits_SNP\n")
  cat(str(subjectHits_SNP))
  cat("\n")
  
  m <- findOverlaps(gr_ROI,gr_TF_EXP_DE)
  
  # cat("m\n")
  # cat(str(m))
  # cat("\n")
  
  subjectHits_ROI<-subjectHits(m)
  
  # cat("subjectHits_ROI\n")
  # cat(str(subjectHits_ROI))
  # cat("\n")
  
  df_SNP <- data.frame(query_region=NA,
                       chr=seqnames(gr_SNP),
                       start=start(gr_SNP),
                       end=end(gr_SNP),
                       strand=c(rep('*', length(gr_SNP))),
                       feature=c(rep('SNP', length(gr_SNP))),
                       source=gr_SNP$name2,
                       Motif_ID=gr_SNP$name3,
                       score=c(rep(0, length(gr_SNP))),
                       Symbol=names(gr_SNP),
                       Intersect_SNP=c(rep('NA', length(gr_SNP))),
                       Intersect_ROI=c(rep('NA', length(gr_SNP))), stringsAsFactors = F)
  # cat("df_SNP_0\n")
  # cat(str(df_SNP))
  # cat("\n")
  
  df_ROI <- data.frame(query_region=NA,
                       chr=seqnames(gr_ROI),
                       start=start(gr_ROI),
                       end=end(gr_ROI),
                       strand=c(rep('*', length(gr_ROI))),
                       feature=names(gr_ROI),
                       source=gr_ROI$name2,
                       Motif_ID=gr_ROI$name3,
                       score=c(rep(0, length(gr_ROI))),
                       Symbol=names(gr_ROI),
                       Intersect_SNP=c(rep('NA', length(gr_ROI))),
                       Intersect_ROI=c(rep('NA', length(gr_ROI))), stringsAsFactors = F)
  
  # cat("df_ROI_0\n")
  # cat(str(df_ROI))
  # cat("\n")
  
  df <- data.frame(query_region=as.character(gr_TF_EXP_DE$name5),
                   chr=seqnames(gr_TF_EXP_DE),
                   start=start(gr_TF_EXP_DE),
                   end=end(gr_TF_EXP_DE),
                   strand=strand(gr_TF_EXP_DE),
                   feature=c(rep('TF_Motif', length(gr_TF_EXP_DE))),
                   source=gr_TF_EXP_DE$name2,
                   Motif_ID=gr_TF_EXP_DE$name3,
                   score=as.numeric(gr_TF_EXP_DE$name4),
                   Symbol=names(gr_TF_EXP_DE),
                   Intersect_SNP=c(rep('NO', length(gr_TF_EXP_DE))),
                   Intersect_ROI=c(rep('NO', length(gr_TF_EXP_DE))), stringsAsFactors = F)
  
  df$Intersect_SNP[subjectHits_SNP]<-'YES'
  df$Intersect_ROI[subjectHits_ROI]<-'YES'
  
  cat("df_0\n")
  cat(str(df))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_ROI))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_ROI)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_SNP))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_SNP)))))
  cat("\n")
  
  check<-df[which(df$Intersect_SNP == 'YES'),]
  
  cat("check_EXP_DE\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(unique(check$Symbol))))
  cat("\n")
  
  
  df<-rbind(df_SNP,df_ROI,df)
  
  df$analysis<-"EXP_DE"
  
  cat("df_1\n")
  cat(str(df))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_ROI))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_ROI)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(df$Intersect_SNP))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(df$Intersect_SNP)))))
  cat("\n")
  
  df_EXP_DE<-df
  
  #### Put together and factorize ----
  
  REP<-unique(rbind(df_EXP,df_EXP_DE))
  
  cat("REP_0\n")
  cat(str(REP))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(REP$source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(REP$source)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(REP$analysis))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(REP$analysis)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(REP$query_region))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(REP$query_region)))))
  cat("\n")
  
  REP$feature<-factor(REP$feature,
                     levels=c("TF_Motif","Del_16bp","Del_80bp","SNP"),
                     ordered=T)

  REP$source<-factor(REP$source,
                     levels=c("CIS_BP","JASPAR","HOMER","ROI","SNP"),
                     ordered=T)
  
  REP$analysis<-factor(REP$analysis,
                     levels=c("EXP","EXP_DE"),
                     ordered=T)
  
  REP$query_region<-factor(REP$query_region,
                       levels=c('7 dna:chromosome chromosome:GRCh37:7:101499894:101499973:1','SNP_change'),
                       ordered=T)
  
  cat("REP_1\n")
  cat(str(REP))
  cat("\n")
  
  
  
  ##### save -----
  
  setwd(out)
 
  write.table(REP, file="TF_motifs_FINAL.tsv", sep="\t", quote=F, row.names = F)
  saveRDS(REP,file="TF_motifs_FINAL.rds")
 
  
  
}

geom_segment_graphs = function(option_list)
{
  
  suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("Biostrings", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("ggseqlogo", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Read input fasta ----
  
  fastaFile<-readDNAStringSet(file = opt$input_fasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df_fasta <- data.frame(seq_name, sequence)
    
 
  cat("df_fasta_0\n")
  cat(str(df_fasta))
  cat("\n")
  
  #### read REP ----
  
  setwd(out)
  
  REP<-readRDS(file="TF_motifs_FINAL.rds")
  REP$Symbol2<-""
  x_SNP<-unique(REP$start[which(REP$feature == 'SNP')])
  
  
  
  cat("REP_0\n")
  cat(str(REP))
  cat("\n")
  cat(sprintf(as.character(names(summary(REP$source)))))
  cat("\n")
  cat(sprintf(as.character(summary(REP$source))))
  cat("\n")
  cat(sprintf(as.character(names(summary(REP$analysis)))))
  cat("\n")
  cat(sprintf(as.character(summary(REP$analysis))))
  cat("\n")
  cat("x_SNP\n")
  cat(sprintf(as.character(x_SNP)))
  cat("\n")
  
  REP<-droplevels(REP[-which(REP$source == 'SNP'),])
  
  cat("REP_1\n")
  cat(str(REP))
  cat("\n")
  
  #### Invert polarity based on strand ----
  
  indx.negative_strand<-which(REP$strand == '-')
  
  cat("indx.negative_strand_1\n")
  cat(str(indx.negative_strand))
  cat("\n")
  
  check<-REP[indx.negative_strand,]
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  NEW_ENDPOINTS<-REP$start[indx.negative_strand]
  NEW_STARTING_POINTS<-REP$end[indx.negative_strand]
  
  REP$start[indx.negative_strand]<-NEW_STARTING_POINTS
  REP$end[indx.negative_strand]<-NEW_ENDPOINTS
  
  check<-REP[indx.negative_strand,]
  
  cat("check_1\n")
  cat(str(check))
  cat("\n")
  
  #### geom segment graphs----
  
  DEBUG <- 1
  
  A_score<-summary(REP$score[!is.na(REP$score)])
  
  
  if(DEBUG == 1)
  {
    cat("A_score\n")
    cat(sprintf(as.character(names(A_score))))
    cat("\n")
    cat(sprintf(as.character(A_score)))
    cat("\n")
  }
  
  max_value<-A_score[6]
  min_value<-A_score[1]
  
  
  step<-round(abs(max_value-min_value)/10,0)
  
  if(step == 0)
  {
    
    step<-1
  }
  breaks.score<-sort(unique(c(0,seq(min_value,max_value+step, by=step))))
  labels.score<-as.character(round(breaks.score),1)
  
  if(DEBUG == 1)
  {
    cat("step_score\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    cat("breaks.score\n")
    cat(sprintf(as.character(breaks.score)))
    cat("\n")
    cat("labels.score\n")
    cat(sprintf(as.character(labels.score)))
    cat("\n")
  }
  
  absolute_min<-min(REP$start[!is.na(REP$start)])
  absolute_max<-max(REP$end[!is.na(REP$end)])
  
  if(DEBUG == 1)
  {
    cat("absolute_min\n")
    cat(sprintf(as.character(absolute_min)))
    cat("\n")
    
    cat("absolute_max\n")
    cat(sprintf(as.character(absolute_max)))
    cat("\n")
    
    
  }
  
  step<-round(abs(absolute_max-absolute_min)/10,0)
  
  if(step == 0)
  {
    
    step<-10
  }
  breaks.COORD<-sort(unique(c(seq(absolute_min,absolute_max+step, by=step))))
  labels.COORD<-as.character(round(breaks.COORD),0)
  
  if(DEBUG == 1)
  {
    cat("step_COORD\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    cat("breaks.COORD\n")
    cat(sprintf(as.character(breaks.COORD)))
    cat("\n")
    cat("labels.COORD\n")
    cat(sprintf(as.character(labels.COORD)))
    cat("\n")
  }
  
  absolute_min<-min(as.numeric(REP$source))
  absolute_max<-max(as.numeric(REP$source))
  
  
  if(DEBUG == 1)
  {
    cat("absolute_min\n")
    cat(sprintf(as.character(absolute_min)))
    cat("\n")
    
    cat("absolute_max\n")
    cat(sprintf(as.character(absolute_max)))
    cat("\n")
   
  }
  
 
  breaks.source<-sort(unique(c(seq(absolute_min,absolute_max, by=1))))
  labels.source<-as.character(c(levels(REP$source)))
  
  
  if(DEBUG == 1)
  {
    cat("step_source\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    cat("breaks.source\n")
    cat(sprintf(as.character(breaks.source)))
    cat("\n")
    cat("labels.source\n")
    cat(sprintf(as.character(labels.source)))
    cat("\n")
  }
  
  if(DEBUG == 1)
  {
    cat("Graph_Part_START:\n")
    
  }
  
  #### Printing LOOP ----
  
  query_region_array<-levels(REP$query_region[!is.na(REP$query_region)])
  
  if(DEBUG == 1)
  {
    cat("query_region_array\n")
    cat(str(query_region_array))
    cat("\n")
  }
  
  for(i in 1:length(query_region_array))
  {
    query_region_array_sel<-query_region_array[i]
    
    cat("------------------query region--------------------------->\t")
    cat(sprintf(as.character(query_region_array_sel)))
    cat("\n")
    
    df_fasta_sel<-df_fasta[which(df_fasta$seq_name == query_region_array_sel),]
    
    if(DEBUG == 1)
    {
      cat("df_fasta_sel_0\n")
      cat(str(df_fasta_sel))
      cat("\n")
      
    }
    
    REP_sel_NA<-unique(REP[is.na(REP$query_region),])
    
    if(DEBUG == 1)
    {
      cat("REP_sel_NA_0\n")
      cat(str(REP_sel_NA))
      cat("\n")
      
    }
    
    REP_sel_region<-unique(REP[which(REP$query_region == query_region_array_sel),])
    
    if(DEBUG == 1)
    {
      cat("REP_sel_region_0\n")
      cat(str(REP_sel_region))
      cat("\n")
      
    }
    
    REP_sel<-rbind(REP_sel_NA,REP_sel_region)
    
    if(DEBUG == 1)
    {
      cat("REP_sel_0\n")
      cat(str(REP_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(REP_sel$source)))))
      cat("\n")
      cat(sprintf(as.character(summary(REP_sel$source))))
      cat("\n")
     
    }
    
    REP_array<-c('HOMER','JASPAR','CIS_BP')
    
    for(k in 1:length(REP_array))
    {
      REP_array_sel<-REP_array[k]
      
      REP_sel$Symbol2<-""
      
      REP_sel$Symbol2[which(REP_sel$source == 'ROI')]<-REP_sel$Symbol[which(REP_sel$source == 'ROI')]
      REP_sel$Symbol2[which(REP_sel$source == REP_array_sel)]<-REP_sel$Symbol[which(REP_sel$source == REP_array_sel)]
      
      
      if(DEBUG == 1)
      {
        check<-REP_sel[which(REP_sel$Symbol == 'PU.1'),]
        
        cat("check_0\n")
        cat(str(check))
        cat("\n")
        cat(sprintf(as.character(check$start)))
        cat("\n")
      }
      
      # jitter_pos_1 <- position_jitter(height = 0.2, seed = 1)
      # 
      # 
      # jitter_pos <- position_jitter(height = 0.4, seed = 1)
      
      jitter_pos_1 <- position_jitter(width = 0, height=0.2, seed = 1)
      
      
      jitter_pos <- position_jitter(width = 0, height=0.4, seed = 1)
      
      
      graph_segment<-ggplot()+
        geom_segment(data=REP_sel[which(REP_sel$source == 'ROI'),],
                     aes(x=start,
                         xend=end,
                         y=as.numeric(source),
                         yend=as.numeric(source)),
                     position=jitter_pos_1,
                     color="black", size=4)+
        geom_segment(data=REP_sel[which(REP_sel$source == 'HOMER'),],
                     aes(x=start,
                         xend=end,
                         y=as.numeric(source),
                         yend=as.numeric(source),
                         color=score),
                     position=jitter_pos,
                     arrow=arrow(angle=40,type="closed",ends="last",length=unit(3,"mm")),
                     size=3)+
        geom_segment(data=REP_sel[which(REP_sel$source == 'JASPAR'),],
                     aes(x=start,
                         xend=end,
                         y=as.numeric(source),
                         yend=as.numeric(source),
                         color=score),
                     position=jitter_pos,
                     arrow=arrow(angle=40,type="closed",ends="last",length=unit(3,"mm")),
                     size=3)+
        geom_segment(data=REP_sel[which(REP_sel$source == 'CIS_BP'),],
                     aes(x=start,
                         xend=end,
                         y=as.numeric(source),
                         yend=as.numeric(source),
                         color=score),
                     position=jitter_pos,
                     arrow=arrow(angle=40,type="closed",ends="last",length=unit(3,"mm")),
                     size=3)+
        scale_color_gradient2(
          low = "blue", 
          mid = "white", 
          high = "red", 
          midpoint = 0,
          breaks=breaks.score,labels=labels.score,
          limits=c(breaks.score[1]-0.01,breaks.score[length(breaks.score)]+0.01),name=paste('logodds',"score",sep="\n"),na.value = "gray")+
        scale_x_continuous(name="chr7", 
                           breaks=breaks.COORD,
                           labels=labels.COORD, 
                           limits=c(breaks.COORD[1],breaks.COORD[length(breaks.COORD)]))+
        scale_y_continuous(name=NULL, 
                           breaks=breaks.source,
                           labels=labels.source, 
                           limits=c(breaks.source[1]-0.5,breaks.source[length(breaks.source)]+0.25))+
        geom_text_repel(data=REP_sel[which(REP_sel$source == 'ROI'),],
                        aes(x=end,
                            y=as.numeric(source),
                            label=Symbol2),
                        position=jitter_pos_1,
                        box.padding = 1,
                        max.overlaps = Inf,
                        show.legend = FALSE)+
        geom_text_repel(data=REP_sel[which(REP_sel$source == 'HOMER'),],
                        aes(x=end,
                            y=as.numeric(source),
                            label=Symbol2),
                        position=jitter_pos,
                        box.padding = 1,
                        max.overlaps = Inf,
                        show.legend = FALSE)+
        geom_text_repel(data=REP_sel[which(REP_sel$source == 'JASPAR'),],
                        aes(x=end,
                            y=as.numeric(source),
                            label=Symbol2),
                        position=jitter_pos,
                        box.padding = 1,
                        max.overlaps = Inf,
                        show.legend = FALSE)+
        geom_text_repel(data=REP_sel[which(REP_sel$source == 'CIS_BP'),],
                        aes(x=end,
                            y=as.numeric(source),
                            label=Symbol2),
                        position=jitter_pos,
                        box.padding = 1,
                        max.overlaps = Inf,
                        show.legend = FALSE)+
        theme_bw()+
        geom_vline(xintercept=x_SNP,linetype="dashed")+
        ggeasy::easy_center_title()
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_MIDDLE:\n")
        
      }
      
      
      
      
      graph_segment<-graph_segment+
        facet_grid(. ~ analysis, scales='free_x', space='free_x') +
        theme_cowplot(font_size = 14)+
        theme( strip.background = element_blank(),
               strip.placement = "inside",
               strip.text = element_text(size=14),
               panel.spacing = unit(0.2, "lines"),
               panel.background=element_rect(fill="white"),
               panel.border=element_rect(colour="black",size=1),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())+
        theme(axis.title.y=element_blank(),
              axis.title.x=element_text(size=24, family="sans"),
              axis.text.y=element_text(angle=90,size=14,color="black", family="sans"),
              axis.text.x=element_text(angle=45, size=12,vjust=1,hjust=1, color="black", family="sans"),
              legend.title=element_text(size=16,color="black", family="sans"),
              legend.text=element_text(size=12,color="black", family="sans"))+
        theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=12))+
        theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
              legend.key.height = unit(1.5, 'cm'), #change legend key height
              legend.key.width = unit(1, 'cm'), #change legend key width
              legend.title = element_text(size=16,color="black", family="sans"), #change legend title font size
              legend.text = element_text(size=14))+ #change legend text font size
        guides(colour = guide_legend(nrow = 1))+
        ggeasy::easy_center_title()
      
      
      
      
      path_graphs<-paste(out,'geom_segment_plots','/',sep='')
      
      if (file.exists(path_graphs)){
        
        
      }else{
        
        dir.create(file.path(path_graphs))
        
      }#path_graphs
      
      
      setwd(path_graphs)
      
      svgname<-paste("Geom_segment_",REP_array_sel,'_',query_region_array_sel,".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= graph_segment,
               device="svg",
               height=10, width=12)
      }
      
      if(DEBUG == 1)
      {
        cat("Graph_Part_END:\n")
        
        
        
      }
      
    }# k in 1:length(REP_array)
    
  }#i in 1:length(query_region_array)
  
  
  
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--input_fasta"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_file_HOMER"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_file_JASPAR"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_file_CIS_BP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_collapsed_EXP_crossed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_collapsed_EXP_and_DE_crossed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ROI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SNP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  # select_motifs(opt)
  geom_segment_graphs(opt)

  
}


###########################################################################

system.time( main() )