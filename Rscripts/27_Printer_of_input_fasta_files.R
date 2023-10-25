
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
suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biostrings", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("seqinr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

data_wrangling = function(option_list)
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
  
  #### READ and transform upstream_span ----
  
  upstream_span = opt$upstream_span
  
  cat("upstream_span_\n")
  cat(sprintf(as.character(upstream_span)))
  cat("\n")
  
  #### READ and transform downstream_span ----
  
  downstream_span = opt$downstream_span
  
  cat("downstream_span_\n")
  cat(sprintf(as.character(downstream_span)))
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
  
  #### READ fastaFile ----
  
  fastaFile<-readDNAStringSet(file = opt$input_fasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df_fasta <- data.frame(seq_name, sequence)
  
  cat("df_fasta_0\n")
  cat(str(df_fasta))
  cat("\n")
  
  df_fasta$seq_name<-gsub("^>","",df_fasta$seq_name)
  df_fasta$chr<-gsub(":.+$","",df_fasta$seq_name)
  df_fasta$start<-gsub("^[^:]+:","",df_fasta$seq_name)
  df_fasta$start<-gsub("-.+$","",df_fasta$start)
  df_fasta$end<-gsub("^[^:]+:[0-9]+-","",df_fasta$seq_name)
  
  cat("df_fasta_1\n")
  cat(str(df_fasta))
  cat("\n")
  
  #### Merge fasta with bed ----
  
  input_bed<-merge(input_bed,
                   df_fasta,
                   by=c('chr','start','end'))
  
  cat("input_bed_POST_merge\n")
  cat(str(input_bed))
  cat("\n")
  
  #### Substring before and after the SNP ----
  
  input_bed$PRE_SNP<-substr(input_bed$sequence, 1, upstream_span)
  input_bed$POST_SNP<-substr(input_bed$sequence, upstream_span+2,nchar(input_bed$sequence))
  
  cat("input_bed_2\n")
  cat(str(input_bed))
  cat("\n")
  
  input_bed$REF_version<-paste(input_bed$PRE_SNP,input_bed$ref,input_bed$POST_SNP, sep='')
  input_bed$ALT_version<-paste(input_bed$PRE_SNP,input_bed$alt,input_bed$POST_SNP, sep='')
  input_bed$name2<-paste(input_bed$rs)
  input_bed$name<-paste(input_bed$seq_name,input_bed$name,sep='|')
  
  cat("input_bed_3\n")
  cat(str(input_bed))
  cat("\n")
  
  
  indx.check<-which(input_bed$REF_version != input_bed$sequence)
  
  cat("indx.check_0\n")
  cat(str(indx.check))
  cat("\n")
  
  check<-input_bed[indx.check,]
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  
  #### subset and do both versions----
  
  
  subset<-unique(input_bed[,c(which(colnames(input_bed) == 'name2'), which(colnames(input_bed) == 'REF_version'), which(colnames(input_bed) == 'ALT_version'))])
  
  cat("subset_0\n")
  cat(str(subset))
  cat("\n")
  
  subset.m<-melt(subset, id.vars='name2', variable.name='version', value.name='sequence')
  
  cat("subset.m_0\n")
  cat(str(subset.m))
  cat("\n")
  
  subset.m$version<-gsub("_version","",subset.m$version)
  
  cat("subset.m_1\n")
  cat(str(subset.m))
  cat("\n")
  
  
  subset.m$version<-factor(subset.m$version,
                           levels=c("REF","ALT"),
                           ordered=T)
  
  cat("subset.m_2\n")
  cat(str(subset.m))
  cat("\n")
  
  subset.m<-subset.m[order(subset.m$name2,subset.m$version),]
  
  
  subset.m$name2<-paste(subset.m$name2,subset.m$version, sep=' ')
  
  cat("subset.m_3\n")
  cat(str(subset.m))
  cat("\n")
  
  
 #### export bed ----
  
  setwd(out)
  
  write.fasta(as.list(subset.m$sequence), subset.m$name, 'ALL_variants_REF_and_ALT.fasta', open = "w", nbchar = 60, as.string = FALSE)
  
  
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
    make_option(c("--input_bed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_fasta"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--upstream_span"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--downstream_span"), type="numeric", default=NULL, 
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
  
  data_wrangling(opt)
  
  
}


###########################################################################

system.time( main() )