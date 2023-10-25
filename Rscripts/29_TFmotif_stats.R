
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

chisq_function = function(option_list)
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
  
  #### READ and transform context_initiator_TFs_subset ----
  
  context_initiator_TFs_subset = unlist(strsplit(opt$context_initiator_TFs_subset, split=","))
  
  cat("context_initiator_TFs_subset_0\n")
  cat(sprintf(as.character(context_initiator_TFs_subset)))
  cat("\n")
  
  #### READ and transform context_only_TFs_subset ----
  
  context_only_TFs_subset = unlist(strsplit(opt$context_only_TFs_subset, split=","))
  
  cat("context_only_TFs_subset_0\n")
  cat(sprintf(as.character(context_only_TFs_subset)))
  cat("\n")
  
  #### READ and transform version_string ----
  
  version_string = unlist(strsplit(opt$version_string, split=","))
  
  cat("version_string_0\n")
  cat(sprintf(as.character(version_string)))
  cat("\n")
  
  #### READ and transform regulation_string ----
  
  regulation_string = unlist(strsplit(opt$regulation_string, split=","))
  
  cat("regulation_string_0\n")
  cat(sprintf(as.character(regulation_string)))
  cat("\n")
  
  regulation_string<-gsub('Transcriptional_R_ATU','Transcriptional R. + ATU',regulation_string)
  regulation_string<-gsub('Transcriptional_R','Transcriptional R.',regulation_string)
  
  cat("regulation_string_1\n")
  cat(sprintf(as.character(regulation_string)))
  cat("\n")
  
  #### READ and transform Intersect_SNP ----
  
  Intersect_SNP = opt$Intersect_SNP
  
  cat("Intersect_SNP_0\n")
  cat(sprintf(as.character(Intersect_SNP)))
  cat("\n")
  
  #### READ Table_S6 ----
  
  Table_S6<-as.data.frame(fread(file=opt$Table_S6, sep="\t", header=T), stringsAsFactors=F)
  
  
  cat("Table_S6_0\n")
  cat(str(Table_S6))
  cat("\n")
  cat(str(unique(Table_S6$VAR_38)))
  cat("\n")
 
  Table_S6_subset<-unique(Table_S6[,c(which(colnames(Table_S6) == 'rs'),
                               which(colnames(Table_S6) == 'VAR'),
                               which(colnames(Table_S6) == 'VAR_38'),
                               which(colnames(Table_S6) == 'Mechanistic_Class'),
                               which(colnames(Table_S6) == 'Manual_curation'),
                               which(colnames(Table_S6) == 'MPRA_CLASS'))])
  
  cat("Table_S6_subset_0\n")
  cat(str(Table_S6_subset))
  cat("\n")
  cat(str(unique(Table_S6_subset$VAR_38)))
  cat("\n")
  
  
  #### READ TF_motifs_annotated ----
  
  TF_motifs_annotated<-as.data.frame(fread(file=opt$TF_motifs_annotated, sep="\t", header=T), stringsAsFactors=F)
  
  
  cat("TF_motifs_annotated_0\n")
  cat(str(TF_motifs_annotated))
  cat("\n")
  cat(str(unique(TF_motifs_annotated$VAR_38)))
  cat("\n")
  
  
  #### merge with annotation ----
  
  
  Table_S6_subset<-merge(Table_S6_subset,
                         TF_motifs_annotated,
                         by="rs")
  
  cat("Table_S6_subset_1\n")
  cat(str(Table_S6_subset))
  cat("\n")
  cat(str(unique(Table_S6_subset$VAR_38)))
  cat("\n")
  
  
  ###### comparison ----
  
  
  df_comparison<-unique(Table_S6_subset[which(Table_S6_subset$version%in%version_string &
                                         Table_S6_subset$Mechanistic_Class%in%regulation_string &
                                         Table_S6_subset$Intersect_SNP%in%Intersect_SNP),])
  
  
  cat("df_comparison_0\n")
  cat(str(df_comparison))
  cat("\n")
  cat(str(unique(df_comparison$VAR_38)))
  cat("\n")
  
  
 
  
  
  indx.context_initiator<-grep(paste(context_initiator_TFs_subset,collapse="|"),df_comparison$Motif_ID)
  
  indx.context_only<-grep(paste(context_only_TFs_subset,collapse="|"),df_comparison$Motif_ID)
  
  
  df_context_initiator<-df_comparison[indx.context_initiator,]
  
  cat("df_context_initiator_0\n")
  cat(str(df_context_initiator))
  cat("\n")
  cat(str(unique(df_context_initiator$VAR_38)))
  cat("\n")
  
  df_context_only<-df_comparison[indx.context_only,]
  
  cat("df_context_only_0\n")
  cat(str(df_context_only))
  cat("\n")
  cat(str(unique(df_context_only$VAR_38)))
  cat("\n")
  
  
  indx.ambivalent<-which(df_context_initiator$rs%in%df_context_only$rs)
  
  cat("indx.ambivalent_0\n")
  cat(str(indx.ambivalent))
  cat("\n")
  
  cat(sprintf(as.character(unique(df_context_initiator$rs[indx.ambivalent]))))
  cat("\n")
  
  df_context_initiator<-df_context_initiator[-indx.ambivalent,]
  
  cat("df_context_initiator_1\n")
  cat(str(df_context_initiator))
  cat("\n")
  cat(str(unique(df_context_initiator$VAR_38)))
  cat("\n")
  
  
  #### Freq_table ----
  
  indx.int<-c(which(colnames(df_context_initiator) =='rs'),which(colnames(df_context_initiator) =='VAR'),which(colnames(df_context_initiator) =='VAR_38'),
              which(colnames(df_context_initiator) =='Mechanistic_Class'),which(colnames(df_context_initiator) =='Manual_curation'),which(colnames(df_context_initiator) =='MPRA_CLASS'))
  
  df_context_initiator_subset<-unique(df_context_initiator[,indx.int])
  
  df_context_initiator_subset$TF_class<-'context_initiator'
  
  cat("df_context_initiator_subset_0\n")
  cat(str(df_context_initiator_subset))
  cat("\n")
  cat(str(unique(df_context_initiator_subset$VAR_38)))
  cat("\n")
  
  indx.int<-c(which(colnames(df_context_only) =='rs'),which(colnames(df_context_only) =='VAR'),which(colnames(df_context_only) =='VAR_38'),
              which(colnames(df_context_only) =='Mechanistic_Class'),which(colnames(df_context_only) =='Manual_curation'),which(colnames(df_context_only) =='MPRA_CLASS'))
  
  df_context_only_subset<-unique(df_context_only[,indx.int])
  
  df_context_only_subset$TF_class<-'context_only'
  
  
  cat("df_context_only_subset_0\n")
  cat(str(df_context_only_subset))
  cat("\n")
  cat(str(unique(df_context_only_subset$VAR_38)))
  cat("\n")
  
  
  Freq_table<-rbind(df_context_initiator_subset,df_context_only_subset)
  
  cat("Freq_table_0\n")
  cat(str(Freq_table))
  cat("\n")
  
  Freq_table.dt<-data.table(Freq_table, key=c("TF_class","MPRA_CLASS"))
  
  
  Freq_table_N<-as.data.frame(Freq_table.dt[,.(Freq=.N), by=key(Freq_table.dt)])
  
  cat("Freq_table_N_0\n")
  cat(str(Freq_table_N))
  cat("\n")
  
  Freq_table_N$MPRA_CLASS<-factor(Freq_table_N$MPRA_CLASS,
                                  levels=c('no MPRA hit','MPRA hit'), ordered=T)
  
  
  Freq_table_N$TF_class<-factor(Freq_table_N$TF_class,
                                  levels=c('context_only','context_initiator'), ordered=T)
  
  cat("Freq_table_N_1\n")
  cat(str(Freq_table_N))
  cat("\n")
  
  Freq_table_N_wide<-as.data.frame(pivot_wider(Freq_table_N,
                                               id_cols=TF_class,
                                               names_from=MPRA_CLASS,
                                               values_from=Freq), stringsAsFactors=F)
  
  cat("Freq_table_N_wide_0\n")
  cat(str(Freq_table_N_wide))
  cat("\n")
  
  Freq_table_N_wide.matrix<-as.matrix(Freq_table_N_wide[,-which(colnames(Freq_table_N_wide) == 'TF_class')])
  
  cat("Freq_table_N_wide.matrix_0\n")
  cat(str(Freq_table_N_wide.matrix))
  cat("\n")
  
  tab.chisq.test<-chisq.test(Freq_table_N_wide.matrix,correct = TRUE)
  
  # cat("tab.chisq.test\n")
  # cat(str(tab.chisq.test))
  # cat("\n")
  pval<-as.numeric(tab.chisq.test$p.value)
  log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
  
  cat("log_pval\n")
  cat(str(log_pval))
  cat("\n")
  
  
  #### SAVE ----
 
  setwd(out)
  
  write.table(file='comparison_table.tsv', df_comparison, sep="\t", quote = F, row.names = F)
  
  
  write.table(rbind(df_context_initiator, df_context_only), file="Selected_variants.tsv", sep = "\t", row.names = F, quote=F)
  
  
  
  Freq_table_N$pval<-pval
  Freq_table_N$log_pval<-log_pval
  
  
  saveRDS(Freq_table_N, file="Chisq_result.rds")
  
}

graph_function = function(option_list)
{
  suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
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
  
  #### READ and transform context_initiator_TFs_subset ----
  
  context_initiator_TFs_subset = unlist(strsplit(opt$context_initiator_TFs_subset, split=","))
  
  cat("context_initiator_TFs_subset_0\n")
  cat(sprintf(as.character(context_initiator_TFs_subset)))
  cat("\n")
  
  #### READ and transform context_only_TFs_subset ----
  
  context_only_TFs_subset = unlist(strsplit(opt$context_only_TFs_subset, split=","))
  
  cat("context_only_TFs_subset_0\n")
  cat(sprintf(as.character(context_only_TFs_subset)))
  cat("\n")
  
  
  #### Read from previous function ----
  
  setwd(out)
  
  filename='Selected_variants.tsv'
  
  Selected_variants<-as.data.frame(fread(file=filename, sep="\t", header=T), stringsAsFactors=F)
  
  cat("Selected_variants_0\n")
  cat(str(Selected_variants))
  cat("\n")
  
  Freq_table_N<-readRDS(file="Chisq_result.rds")
  
  cat("Freq_table_N_0\n")
  cat(str(Freq_table_N))
  cat("\n")
  
  Freq_table_N$TFs<-NA
  
  
  Freq_table_N$TFs[which(Freq_table_N$TF_class =='context_only')]<-paste(context_only_TFs_subset, collapse=',')
  Freq_table_N$TFs[which(Freq_table_N$TF_class =='context_initiator')]<-paste(context_initiator_TFs_subset, collapse=',')
  
  
  cat("Freq_table_N_1\n")
  cat(str(Freq_table_N))
  cat("\n")
  
  #### graph ----
  
  
  breaks.Rank<-(seq(0,max(Freq_table_N$Freq)+2,by=1))
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  vector_colors<-brewer.pal(8,"Paired")[c(8,7)]
  
  cat("vector_colors_0\n")
  cat(str(vector_colors))
  cat("\n")
  
  
  stacked_barplot<-ggplot(data=Freq_table_N,
                          aes(x=MPRA_CLASS, y=Freq, fill=TF_class, label=TFs)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=16, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"))+
    theme(axis.text.x=element_text(angle=45,size=14,vjust=1,hjust=1,color="black", family="sans"))+
    scale_y_continuous(name=paste('Variants with Transcriptional R. or Transcriptional R. + ATU',sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+0.1))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_fill_manual(values=vector_colors,drop=F)+
    theme(legend.position="right",legend.title=element_blank(), legend.text = element_text(size=14))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    geom_text_repel(aes(color=TF_class),
                    position = position_nudge_repel(x = 0.2, y = -0.2),
                    box.padding = 1,
                    size=4)+
    scale_color_manual(values=vector_colors,drop=F)+
    ggeasy::easy_center_title()
  
  
  
  setwd(out)
  
  svgname<-paste(paste("Stacked_barplot_TF_motifs", sep='_'),".svg",sep='')
  
  
  svglite(svgname, width = 8, height = 8)
  print(stacked_barplot)
  dev.off()
  
  
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
    make_option(c("--Table_S6"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_motifs_annotated"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Intersect_SNP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--version_string"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--regulation_string"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--context_only_TFs_subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--context_initiator_TFs_subset"), type="character", default=NULL, 
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
  
  chisq_function(opt)
  graph_function(opt)
  
  
}


###########################################################################

system.time( main() )