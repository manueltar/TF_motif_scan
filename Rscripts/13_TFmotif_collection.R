
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
                                which(colnames(input_file_HOMER) == 'source'))])
  
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
                                       which(colnames(input_file_JASPAR_DEF) == 'source'))])
  
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
                                       which(colnames(input_file_CIS_BP) == 'source'))])
  
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
  
  ALL_df.dt<-data.table(ALL_df, key=c("Symbol","query_region"))
  
  
  ALL_df_collapsed<-as.data.frame(ALL_df.dt[,.(string_source=paste(source, collapse = ";")), by=key(ALL_df.dt)], stringsAsFactors=F)
  
  cat("ALL_df_collapsed_0\n")
  cat(str(ALL_df_collapsed))
  cat("\n")
 
  ALL_df_collapsed$ensembl_gene_id <- mapIds(org.Hs.eg.db, keys=ALL_df_collapsed$Symbol, keytype="SYMBOL",
                           column="ENSEMBL", multiVals="first")
  
  cat("ALL_df_collapsed_1\n")
  cat(str(ALL_df_collapsed))
  cat("\n")
  
  ### harcoded rcovery of PU.1 and others ----
  
  ALL_df_collapsed$ensembl_gene_id[which(ALL_df_collapsed$Symbol == 'PU.1')]<-'ENSG00000066336'
  ALL_df_collapsed$ensembl_gene_id[which(ALL_df_collapsed$Symbol == 'E2A')]<-'ENSG00000071564'
  ALL_df_collapsed$ensembl_gene_id[which(ALL_df_collapsed$Symbol == 'HEB')]<-'ENSG00000140262'
  ALL_df_collapsed$ensembl_gene_id[which(ALL_df_collapsed$Symbol == 'RARa')]<-'ENSG00000131759'
  ALL_df_collapsed$ensembl_gene_id[which(ALL_df_collapsed$Symbol == 'TCFL2')]<-'ENSG00000099949'
  ALL_df_collapsed$ensembl_gene_id[which(ALL_df_collapsed$Symbol == 'THRb')]<-'ENSG00000151090'
  
  
  #### NO NA collapse ----
  
  ALL_df_collapsed_NO_NA<-ALL_df_collapsed[!is.na(ALL_df_collapsed$ensembl_gene_id),]
  
  cat("ALL_df_collapsed_NO_NA_0\n")
  cat(str(ALL_df_collapsed_NO_NA))
  cat("\n")
  
  
  
  ##### save -----
  setwd(out)
  write.table(ALL_df_collapsed_NO_NA, file="TF_motifs_collapsed.tsv", sep="\t", quote=F, row.names = F)
  
  # cat("Hello_world\n")
  # write.table(ALL_df_collapsed, file="test.tsv", sep="\t", quote=F, row.names = F)
  
  
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
    make_option(c("--input_file_HOMER"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_file_JASPAR"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_file_CIS_BP"), type="character", default=NULL, 
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

  
}


###########################################################################

system.time( main() )