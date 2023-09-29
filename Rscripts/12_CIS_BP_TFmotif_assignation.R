
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


opt = NULL

options(warn = 1)

Assign_TF = function(option_list)
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
  
  #### READ and transform output_CIS_BP_assigned ----
  
  output_CIS_BP_assigned = opt$output_CIS_BP_assigned
  
  cat("output_CIS_BP_assigned_\n")
  cat(sprintf(as.character(output_CIS_BP_assigned)))
  cat("\n")
  
  
  #### READ and transform TF_Species ----
  
  TF_Species = opt$TF_Species
  
  cat("TF_Species_\n")
  cat(sprintf(as.character(TF_Species)))
  cat("\n")
  
  #### READ and transform TF_Status ----
  
  TF_Status = opt$TF_Status
  
  cat("TF_Status_\n")
  cat(sprintf(as.character(TF_Status)))
  cat("\n")
  
  #### READ input_file ----
  
  TF_file<-as.data.frame(fread(file=opt$input_file, sep="\t", header=F, skip=5), stringsAsFactors=F)
  colnames(TF_file)<-c('query_region','start','end','Motif_ID','score','strand')
  
  cat("TF_file_0\n")
  cat(str(TF_file))
  cat("\n")
  cat(str(unique(TF_file$Motif_ID)))
  cat("\n")
  
  #### READ TF_key_library ----
  
  TF_key_library<-as.data.frame(fread(file=opt$TF_key_library, sep="\t", header=T), stringsAsFactors=F)
  
  cat("TF_key_library_0\n")
  cat(str(TF_key_library))
  cat("\n")
  
  TF_key_library_filtered<-TF_key_library[which(TF_key_library$TF_Species == TF_Species &
                                                  TF_key_library$TF_Status ==  TF_Status),]
  
  cat("TF_key_library_filtered_0\n")
  cat(str(TF_key_library_filtered))
  cat("\n")
  cat(str(unique(TF_key_library_filtered$Motif_ID)))
  cat("\n")
  cat(str(unique(TF_key_library_filtered$DBID)))
  cat("\n")
  
  #### Merge ----
  
  
  TF_file<-merge(TF_file,
                 TF_key_library_filtered,
                 by='Motif_ID')
  
  cat("TF_file_0\n")
  cat(str(TF_file))
  cat("\n")
  cat(str(unique(TF_file$Motif_ID)))
  cat("\n")
  cat(str(unique(TF_file$DBID)))
  cat("\n")
  
  write.table(TF_file,file=output_CIS_BP_assigned, sep="\t", quote=F,row.names = F)
 
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
    make_option(c("--input_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_key_library"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_Species"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_Status"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--output_CIS_BP_assigned"), type="character", default=NULL, 
                metavar="filename", 
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
  
  Assign_TF(opt)

  
}


###########################################################################

system.time( main() )