
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

Cross_function = function(option_list)
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
  
 
  #### READ TF_collapsed ----
  
  TF_collapsed<-as.data.frame(fread(file=opt$TF_collapsed, sep="\t", header=T), stringsAsFactors=F)

  cat("TF_collapsed_0\n")
  cat(str(TF_collapsed))
  cat("\n")
  cat(str(unique(TF_collapsed$Symbol)))
  cat("\n")
  cat(str(unique(TF_collapsed$ensembl_gene_id)))
  cat("\n")
  
  #### Read post normalised & filtered matrix----
  
  POST_NOR_subset<-as.data.frame(fread(file=opt$POST_NOR_subset, sep="\t"), stringsAsFactors=F)
  colnames(POST_NOR_subset)[which(colnames(POST_NOR_subset) == 'V1')]<-'ensembl_gene_id'
  
  cat("POST_NOR_subset_0\n")
  cat(str(POST_NOR_subset))
  cat("\n")
  
  #### Read DE genes ----
  
  DE_edgeR_results<-as.data.frame(fread(file=opt$DE_edgeR_results, sep="\t"), stringsAsFactors=F)

  cat("DE_edgeR_results_0\n")
  cat(str(DE_edgeR_results))
  cat("\n")
  
  DE_edgeR_results_ctrl<-unique(DE_edgeR_results[which(DE_edgeR_results$comparison == 'Time_course_wt'),])
  
  
  cat("DE_edgeR_results_ctrl_0\n")
  cat(str(DE_edgeR_results_ctrl))
  cat("\n")
  
  DE_edgeR_results_ctrl_SIG<-DE_edgeR_results_ctrl[which(DE_edgeR_results_ctrl$Minus_logpval >=1.3),]
  
  cat("DE_edgeR_results_ctrl_SIG_0\n")
  cat(str(DE_edgeR_results_ctrl_SIG))
  cat("\n")
  
  #### select TFs expressed in our cells ----
  
  TF_collapsed_subset<-unique(TF_collapsed[which(TF_collapsed$ensembl_gene_id%in%POST_NOR_subset$ensembl_gene_id),])
  
  cat("TF_collapsed_subset_0\n")
  cat(str(TF_collapsed_subset))
  cat("\n")
  cat(str(unique(TF_collapsed_subset$Symbol)))
  cat("\n")
  cat(str(unique(TF_collapsed_subset$ensembl_gene_id)))
  cat("\n")
  
  TF_collapsed_subset_DE<-TF_collapsed_subset[which(TF_collapsed_subset$ensembl_gene_id%in%DE_edgeR_results_ctrl_SIG$ensembl_gene_id),]
  
  cat("TF_collapsed_subset_DE_0\n")
  cat(str(TF_collapsed_subset_DE))
  cat("\n")
  cat(str(unique(TF_collapsed_subset_DE$Symbol)))
  cat("\n")
  cat(str(unique(TF_collapsed_subset_DE$ensembl_gene_id)))
  cat("\n")
  
  ##### save -----
  setwd(out)
  write.table(TF_collapsed_subset, file="TF_motifs_collapsed_crossed_with_K562_EXP.tsv", sep="\t", quote=F, row.names = F)
  write.table(TF_collapsed_subset_DE, file="TF_motifs_collapsed_crossed_with_K562_EXP_DE_in_ctrl.tsv", sep="\t", quote=F, row.names = F)
  
  
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
    make_option(c("--TF_collapsed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--POST_NOR_subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DE_edgeR_results"), type="character", default=NULL, 
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
  
  Cross_function(opt)

  
}


###########################################################################

system.time( main() )