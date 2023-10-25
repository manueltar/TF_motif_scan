
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
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

collect = function(option_list)
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
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ input_bed ----
  
  input_bed<-as.data.frame(fread(file=opt$input_bed, sep="\t", header=F), stringsAsFactors=F)
  colnames(input_bed)<-c('chr','start','end','name','score','strand')
  
  cat("input_bed_0\n")
  cat(str(input_bed))
  cat("\n")
  
  input_bed$rs<-gsub("__.+$","",input_bed$name)
  input_bed$VAR_38<-gsub("^[^__]+__","",input_bed$name)
  input_bed$ref<-gsub("^[^_]+_[^_]+_","",input_bed$VAR_38)
  input_bed$ref<-gsub("_.+$","",input_bed$ref)
  input_bed$alt<-gsub("^[^_]+_[^_]+_[^_]+_","",input_bed$VAR_38)
  
  cat("input_bed_1\n")
  cat(str(input_bed))
  cat("\n")
 
  
  ind.dep<-c(which(colnames(input_bed) == 'name'),which(colnames(input_bed) == 'strand'),which(colnames(input_bed) == 'score'))
  
  
  input_bed_subset<-unique(input_bed[,-ind.dep])
  
  colnames(input_bed_subset)[which(colnames(input_bed_subset) == 'start')]<-'Region_start'
  colnames(input_bed_subset)[which(colnames(input_bed_subset) == 'end')]<-'Region_end'
  
  
  cat("input_bed_subset_0\n")
  cat(str(input_bed_subset))
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
  
  df_HOMER<-unique(input_file_HOMER[,c(which(colnames(input_file_HOMER) == 'query_region'),
                                which(colnames(input_file_HOMER) == 'Symbol'),
                                which(colnames(input_file_HOMER) == 'source'),
                                which(colnames(input_file_HOMER) == 'start'),
                                which(colnames(input_file_HOMER) == 'end'),
                                which(colnames(input_file_HOMER) == 'Motif_ID'),
                                which(colnames(input_file_HOMER) == 'strand'),
                                which(colnames(input_file_HOMER) == 'score'))])
  
  cat("df_HOMER_0\n")
  cat(str(df_HOMER))
  cat("\n")
  cat(str(unique(df_HOMER$Symbol)))
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
  
  df_JASPAR<-unique(input_file_JASPAR_DEF[,c(which(colnames(input_file_JASPAR_DEF) == 'query_region'),
                                             which(colnames(input_file_JASPAR_DEF) == 'Symbol'),
                                             which(colnames(input_file_JASPAR_DEF) == 'source'),
                                             which(colnames(input_file_JASPAR_DEF) == 'start'),
                                             which(colnames(input_file_JASPAR_DEF) == 'end'),
                                             which(colnames(input_file_JASPAR_DEF) == 'Motif_ID'),
                                             which(colnames(input_file_JASPAR_DEF) == 'strand'),
                                             which(colnames(input_file_JASPAR_DEF) == 'score'))])
  
  cat("df_JASPAR_0\n")
  cat(str(df_JASPAR))
  cat("\n")
  cat(str(unique(df_JASPAR$Symbol)))
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
  
  df_CIS_BP<-unique(input_file_CIS_BP[,c(which(colnames(input_file_CIS_BP) == 'query_region'),
                                         which(colnames(input_file_CIS_BP) == 'Symbol'),
                                         which(colnames(input_file_CIS_BP) == 'source'),
                                         which(colnames(input_file_CIS_BP) == 'start'),
                                         which(colnames(input_file_CIS_BP) == 'end'),
                                         which(colnames(input_file_CIS_BP) == 'Motif_ID'),
                                         which(colnames(input_file_CIS_BP) == 'strand'),
                                         which(colnames(input_file_CIS_BP) == 'score'))])
  
  cat("df_CIS_BP_0\n")
  cat(str(df_CIS_BP))
  cat("\n")
  cat(str(unique(df_CIS_BP$Symbol)))
  cat("\n")
 
  
  #### Merge all ----
  
  ALL_df<-rbind(df_HOMER,df_JASPAR,df_CIS_BP)
  
  cat("ALL_df_0\n")
  cat(str(ALL_df))
  cat("\n")
  cat(str(unique(ALL_df$Symbol)))
  cat("\n")
  
  
  ALL_df$ensembl_gene_id <- mapIds(org.Hs.eg.db, keys=ALL_df$Symbol, keytype="SYMBOL",
                           column="ENSEMBL", multiVals="first")
  
  cat("ALL_df_1\n")
  cat(str(ALL_df))
  cat("\n")
  
  ### harcoded recovery of PU.1 and others due to alias problems with the motif dBs ----
  
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'PU.1')]<-'ENSG00000066336'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'PU.1-IRF')]<-'ENSG00000066336'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'PU.1:IRF8')]<-'ENSG00000066336'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'E2A')]<-'ENSG00000071564'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'HEB')]<-'ENSG00000140262'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'RARa')]<-'ENSG00000131759'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'TCFL2')]<-'ENSG00000099949'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'THRb')]<-'ENSG00000151090'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'BMYB')]<-'ENSG00000101057'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'HIF-1a')]<-'ENSG00000100644'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'HIF-1b')]<-'ENSG00000143437'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'EAR2')]<-'ENSG00000160113'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'NFkB-p50,p52')]<-'ENSG00000109320,ENSG00000077150'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'RXR')]<-'ENSG00000186350'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'ZNF143|STAF')]<-'ENSG00000166478'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'AP-1')]<-'ENSG00000175592'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'TEAD')]<-'ENSG00000187079'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'EKLF')]<-'ENSG00000105610'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'AMYB')]<-'ENSG00000185697'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'HNF6')]<-'ENSG00000169856'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'ETS')]<-'ENSG00000134954'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol%in%c('NFkB-p65-Rel','NFkB-p65'))]<-'ENSG00000173039'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'HIF2a')]<-'ENSG00000116016'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'NPAS')]<-'ENSG00000130751'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'bHLHE40')]<-'ENSG00000134107'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'bHLHE41')]<-'ENSG00000123095'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'BMAL1')]<-'ENSG00000133794'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'RORa')]<-'ENSG00000069667'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'RORg')]<-'ENSG00000143365'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'NFAT')]<-'ENSG00000131196'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'CEBP')]<-'ENSG00000245848'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'RUNX')]<-'ENSG00000159216'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'LXRE')]<-'ENSG00000025434'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'NFY')]<-'ENSG00000001167'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'STAT5')]<-'ENSG00000126561'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'E2F')]<-'ENSG00000101412'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'RARg')]<-'ENSG00000172819'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'EBF')]<-'ENSG00000164330'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'CEBP:CEBP')]<-'ENSG00000245848'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'BORIS')]<-'ENSG00000124092'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'CENPBD1')]<-'ENSG00000177946'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'ZSCAN5')]<-'ENSG00000131848'
  ALL_df$ensembl_gene_id[which(ALL_df$Symbol == 'HKR1')]<-'ENSG00000181666'
  
  
  
  
  
  
  
  ALL_df$rs<-gsub(" .+$","",ALL_df$query_region)
  ALL_df$version<-gsub("^[^ ]+ ","",ALL_df$query_region)
  
  
  cat("ALL_df_2\n")
  cat(str(ALL_df))
  cat("\n")
  
  
  check<-ALL_df[is.na(ALL_df$ensembl_gene_id),]
  
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$Symbol)))
  cat("\n")
  cat(sprintf(as.character(unique(check$Symbol))))
  cat("\n")
  
  #### bind the input bed subset ----
  
  ALL_df_NO_NA<-ALL_df[!is.na(ALL_df$ensembl_gene_id),]
  
  cat("ALL_df_NO_NA_0\n")
  cat(str(ALL_df_NO_NA))
  cat("\n")
  cat(str(unique(ALL_df_NO_NA$Symbol)))
  cat("\n")
  
  ALL_df_NO_NA<-merge(input_bed_subset,
                          ALL_df_NO_NA,
                          by="rs",
                          all=T)
  
  ALL_df_NO_NA$start<-ALL_df_NO_NA$Region_start+ALL_df_NO_NA$start+1
  ALL_df_NO_NA$end<-ALL_df_NO_NA$Region_start+ALL_df_NO_NA$end
  
  
  cat("ALL_df_NO_NA_1\n")
  cat(str(ALL_df_NO_NA))
  cat("\n")
  cat(str(unique(ALL_df_NO_NA$Symbol)))
  cat("\n")
 
  ##### save -----
  
  setwd(out)
  write.table(ALL_df_NO_NA, file="TF_motifs_collapsed.tsv", sep="\t", quote=F, row.names = F)
  
  # cat("Hello_world\n")
  # write.table(ALL_df, file="test.tsv", sep="\t", quote=F, row.names = F)
  
  
}

intersect_with_SNP = function(option_list)
{
  suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
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
  
  #### READ and transform context_initiator_TFs_string ----
  
  context_initiator_TFs_string = unlist(strsplit(opt$context_initiator_TFs_string, split=","))
  
  cat("context_initiator_TFs_string_0\n")
  cat(sprintf(as.character(context_initiator_TFs_string)))
  cat("\n")
  
  #### READ input from the previous step  ----
  
  setwd(out)
  
  filename<-"TF_motifs_collapsed.tsv"
  
  ALL_df_NO_NA<-as.data.frame(fread(file=filename, sep="\t", header=T))
  
  cat("ALL_df_NO_NA_1\n")
  cat(str(ALL_df_NO_NA))
  cat("\n")
  cat(str(unique(ALL_df_NO_NA$Symbol)))
  cat("\n")
  
  
  VARS_df <- unique(data.frame(chr=gsub("_.+$","",ALL_df_NO_NA$VAR_38),
                        pos38=gsub("^[^_]+_","",ALL_df_NO_NA$VAR_38),
                        names=ALL_df_NO_NA$VAR_38, stringsAsFactors = F))
  
  VARS_df$pos38<-as.integer(gsub("_.+$","",VARS_df$pos38))
  
  cat("VARS_df_0\n")
  cat(str(VARS_df))
  cat("\n")
  
  
  gr_VARS <- GRanges(
    seqnames = as.character(VARS_df$chr),
    strand="*",
    ranges=IRanges(
      start=VARS_df$pos38,
      names=VARS_df$names,
      end=VARS_df$pos38))
  
  cat("gr_VARS_0\n")
  cat(str(gr_VARS))
  cat("\n")
  
  #### HOMER ----
  
  ALL_df_NO_NA_HOMER<-droplevels(ALL_df_NO_NA[which(ALL_df_NO_NA$source == 'HOMER'),])
  
  cat("ALL_df_NO_NA_HOMER_0\n")
  cat(str(ALL_df_NO_NA_HOMER))
  cat("\n")
  cat(str(unique(ALL_df_NO_NA_HOMER$Symbol)))
  cat("\n")
  
  gr_HOMER <- GRanges(
    seqnames = as.character(ALL_df_NO_NA_HOMER$chr),
    name2='HOMER',
    name3=ALL_df_NO_NA_HOMER$Motif_ID,
    name4=ALL_df_NO_NA_HOMER$score,
    name5=ALL_df_NO_NA_HOMER$rs,
    name6=ALL_df_NO_NA_HOMER$version,
    name7=ALL_df_NO_NA_HOMER$ensembl_gene_id,
    strand=ALL_df_NO_NA_HOMER$strand,
    ranges=IRanges(
      start=ALL_df_NO_NA_HOMER$start,
      end=ALL_df_NO_NA_HOMER$end,
      names=ALL_df_NO_NA_HOMER$Symbol))
  
  cat("gr_HOMER_0\n")
  cat(str(gr_HOMER))
  cat("\n")
  
  #### JASPAR ----
  
  ALL_df_NO_NA_JASPAR<-droplevels(ALL_df_NO_NA[which(ALL_df_NO_NA$source == 'JASPAR'),])
  
  cat("ALL_df_NO_NA_JASPAR_0\n")
  cat(str(ALL_df_NO_NA_JASPAR))
  cat("\n")
  cat(str(unique(ALL_df_NO_NA_JASPAR$Symbol)))
  cat("\n")
  
  gr_JASPAR <- GRanges(
    seqnames = as.character(ALL_df_NO_NA_JASPAR$chr),
    name2='JASPAR',
    name3=ALL_df_NO_NA_JASPAR$Motif_ID,
    name4=ALL_df_NO_NA_JASPAR$score,
    name5=ALL_df_NO_NA_JASPAR$rs,
    name6=ALL_df_NO_NA_JASPAR$version,
    name7=ALL_df_NO_NA_JASPAR$ensembl_gene_id,
    strand=ALL_df_NO_NA_JASPAR$strand,
    ranges=IRanges(
      start=ALL_df_NO_NA_JASPAR$start,
      end=ALL_df_NO_NA_JASPAR$end,
      names=ALL_df_NO_NA_JASPAR$Symbol))
  
  cat("gr_JASPAR_0\n")
  cat(str(gr_JASPAR))
  cat("\n")
  
  #### CIS_BP ----
  
  ALL_df_NO_NA_CIS_BP<-droplevels(ALL_df_NO_NA[which(ALL_df_NO_NA$source == 'CIS_BP'),])
  
  cat("ALL_df_NO_NA_CIS_BP_0\n")
  cat(str(ALL_df_NO_NA_CIS_BP))
  cat("\n")
  cat(str(unique(ALL_df_NO_NA_CIS_BP$Symbol)))
  cat("\n")
  
  gr_CIS_BP <- GRanges(
    seqnames = as.character(ALL_df_NO_NA_CIS_BP$chr),
    name2='CIS_BP',
    name3=ALL_df_NO_NA_CIS_BP$Motif_ID,
    name4=ALL_df_NO_NA_CIS_BP$score,
    name5=ALL_df_NO_NA_CIS_BP$rs,
    name6=ALL_df_NO_NA_CIS_BP$version,
    name7=ALL_df_NO_NA_CIS_BP$ensembl_gene_id,
    strand=ALL_df_NO_NA_CIS_BP$strand,
    ranges=IRanges(
      start=ALL_df_NO_NA_CIS_BP$start,
      end=ALL_df_NO_NA_CIS_BP$end,
      names=ALL_df_NO_NA_CIS_BP$Symbol))
  
  cat("gr_CIS_BP_0\n")
  cat(str(gr_CIS_BP))
  cat("\n")
  
  gr_TF<-append(gr_HOMER,gr_JASPAR)
  # gr_TF<-append(gr_TF,gr_CIS_BP) # Dont add gr_CIS_BP
  
  #### find overlap with SNP ---
  
  m <- findOverlaps(gr_VARS,gr_TF)
  
  # cat("m\n")
  # cat(str(m))
  # cat("\n")
  
  subjectHits_SNP<-subjectHits(m)
  
  cat("subjectHits_SNP\n")
  cat(str(subjectHits_SNP))
  cat("\n")
  
  df <- data.frame(rs=as.character(gr_TF$name5),
                   version=as.character(gr_TF$name6),
                   chr=as.character(seqnames(gr_TF)),
                   start=as.integer(start(gr_TF)),
                   end=as.integer(end(gr_TF)),
                   strand=as.character(strand(gr_TF)),
                   source=as.character(gr_TF$name2),
                   Motif_ID=as.character(gr_TF$name3),
                   score=as.numeric(gr_TF$name4),
                   Symbol=names(gr_TF),
                   ensembl_gene_id=as.character(gr_TF$name7),
                   Intersect_SNP=c(rep('NO', length(as.character(seqnames(gr_TF))))), stringsAsFactors = F)
  
  cat("df_0\n")
  cat(str(df))
  cat("\n")
  
  df$Intersect_SNP[subjectHits_SNP]<-'YES'
  
  cat("df_1\n")
  cat(str(df))
  cat("\n")
  
  check<-df[which(df$Intersect_SNP == 'YES' &
                    df$version == 'REF'),]
  
  cat("check\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$Symbol)))
  cat("\n")
  cat(str(unique(check$rs)))
  cat("\n")
  
  #### context_initiator_TFs_string ----
  
  indx.context_initiator<-grep(paste(context_initiator_TFs_string,collapse="|"),df$Motif_ID)
  
  df$context_initiator<-'NO'
  
  df$context_initiator[indx.context_initiator]<-'YES'
  
  check2<-df[which(df$context_initiator == 'YES' &
                     df$version == 'REF'),]
  
  cat("check2\n")
  cat(str(check2))
  cat("\n")
  cat(str(unique(check2$Symbol)))
  cat("\n")
  cat(str(unique(check2$rs)))
  cat("\n")
  cat(sprintf(as.character(unique(check2$Symbol))))
  cat("\n")
  
  #### check3 ----
  
  check3<-df[which(df$context_initiator == 'YES' &
                     df$Intersect_SNP == 'YES' &
                     df$version == 'REF'),]
  
  cat("check3\n")
  cat(str(check3))
  cat("\n")
  cat(str(unique(check3$Symbol)))
  cat("\n")
  cat(str(unique(check3$rs)))
  cat("\n")
  cat(sprintf(as.character(unique(check3$Symbol))))
  cat("\n")
  
  
  ##### save -----
  
  setwd(out)
  
  write.table(df, file="TF_motifs_FINAL.tsv", sep="\t", quote=F, row.names = F)
  saveRDS(df,file="TF_motifs_FINAL.rds")
  
  
  # write.table(check3, file="test.tsv", sep="\t", quote=F, row.names = F)
  
  
}

Add_EXP_data = function(option_list)
{
  suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
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
  
  #### Read post normalised & filtered matrix----
  
  POST_NOR_subset<-as.data.frame(fread(file=opt$POST_NOR_subset, sep="\t"), stringsAsFactors=F)
  colnames(POST_NOR_subset)[which(colnames(POST_NOR_subset) == 'V1')]<-'ensembl_gene_id'
  
  cat("POST_NOR_subset_0\n")
  cat(str(POST_NOR_subset))
  cat("\n")
  
  #### READ input from the previous step  ----
  
  setwd(out)
  
  filename<-"TF_motifs_FINAL.tsv"
  
  TF_motifs_FINAL<-as.data.frame(fread(file=filename, sep="\t", header=T))
  
  cat("TF_motifs_FINAL_1\n")
  cat(str(TF_motifs_FINAL))
  cat("\n")
  cat(str(unique(TF_motifs_FINAL$Symbol)))
  cat("\n")
  
  
  TF_motifs_FINAL$K562_EXP<-NA
  
  TF_motifs_FINAL$K562_EXP[which(TF_motifs_FINAL$ensembl_gene_id%in%POST_NOR_subset$ensembl_gene_id)]<-'YES'
  TF_motifs_FINAL$K562_EXP[-which(TF_motifs_FINAL$ensembl_gene_id%in%POST_NOR_subset$ensembl_gene_id)]<-'NO'
  
  #### check3 ----
  
  check3<-TF_motifs_FINAL[which(TF_motifs_FINAL$context_initiator == 'YES' &
                     TF_motifs_FINAL$Intersect_SNP == 'YES' &
                     TF_motifs_FINAL$version == 'REF' &
                     TF_motifs_FINAL$K562_EXP =='YES'),]
  
  cat("check3\n")
  cat(str(check3))
  cat("\n")
  cat(str(unique(check3$Symbol)))
  cat("\n")
  cat(str(unique(check3$rs)))
  cat("\n")
  cat(sprintf(as.character(unique(check3$Symbol))))
  cat("\n")
  
  #### SAVE ----
  
  setwd(out)
  
  write.table(TF_motifs_FINAL, file="TF_motifs_FINAL_K562_EXP.tsv", sep="\t", quote=F, row.names = F)
  # write.table(check3, file="test.tsv", sep="\t", quote=F, row.names = F)
  
  
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
    make_option(c("--POST_NOR_subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--context_initiator_TFs_string"), type="character", default=NULL, 
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
    make_option(c("--input_bed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
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
  
  collect(opt)
  intersect_with_SNP(opt)
  Add_EXP_data(opt)

  
}


###########################################################################

system.time( main() )