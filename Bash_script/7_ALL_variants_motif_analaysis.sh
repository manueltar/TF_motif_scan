#!/bin/bash

analysis=$1
MASTER_ROUTE=$2
upstream_span=$3
downstream_span=$4

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0



 output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

rm -rf $output_dir
mkdir -p $output_dir

 Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files


#### bedtools_getfasta_intervals #############################


type=$(echo "bedtools_getfasta_intervals""_""$analysis")
outfile_bedtools_getfasta_intervals=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_bedtools_getfasta_intervals
echo -n "" > $outfile_bedtools_getfasta_intervals
name_bedtools_getfasta_intervals=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_bedtools_getfasta_intervals=$(echo "$Rscripts_path""26_printer_bedtools_getfasta_interval_bed.R")


Table_S6=$(echo "/group/soranzo/manuel.tardaguila/Motif_analysis/""Table_S6.tsv")

myjobid_bedtools_getfasta_intervals=$(sbatch --job-name=$name_bedtools_getfasta_intervals --output=$outfile_bedtools_getfasta_intervals --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_bedtools_getfasta_intervals --Table_S6 $Table_S6 --type $type --upstream_span $upstream_span --downstream_span $downstream_span --out $output_dir")
myjobid_seff_bedtools_getfasta_intervals=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta_intervals --open-mode=append --output=$outfile_bedtools_getfasta_intervals --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bedtools_getfasta_intervals >> $outfile_bedtools_getfasta_intervals")


#### bedtools_getfasta #############################


type=$(echo "bedtools_getfasta""_""$analysis")
outfile_bedtools_getfasta=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_bedtools_getfasta
echo -n "" > $outfile_bedtools_getfasta
name_bedtools_getfasta=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

module load bedtools2/2.31.0

reference_genome=$(echo "/processing_data/reference_datasets/iGenomes/2022.1/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa")
input_bed=$(echo "$output_dir""VARS.bed")
output_fasta=$(echo "$output_dir""intervals.fasta")

myjobid_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta_intervals --job-name=$name_bedtools_getfasta --output=$outfile_bedtools_getfasta --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="bedtools getfasta -fi $reference_genome -bed $input_bed > $output_fasta")
myjobid_seff_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --open-mode=append --output=$outfile_bedtools_getfasta --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bedtools_getfasta >> $outfile_bedtools_getfasta")


#### printer_fasta_files #############################


type=$(echo "printer_fasta_files""_""$analysis")
outfile_printer_fasta_files=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_printer_fasta_files
echo -n "" > $outfile_printer_fasta_files
name_printer_fasta_files=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_printer_fasta_files=$(echo "$Rscripts_path""27_Printer_of_input_fasta_files.R")


input_bed=$(echo "$output_dir""VARS.bed")
input_fasta=$(echo "$output_dir""intervals.fasta")


myjobid_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --job-name=$name_printer_fasta_files --output=$outfile_printer_fasta_files --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_printer_fasta_files --input_bed $input_bed --input_fasta $input_fasta --upstream_span $upstream_span --downstream_span $downstream_span --type $type --out $output_dir")
myjobid_seff_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --open-mode=append --output=$outfile_printer_fasta_files --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_printer_fasta_files >> $outfile_printer_fasta_files")

#myjobid_printer_fasta_files=$(sbatch --job-name=$name_printer_fasta_files --output=$outfile_printer_fasta_files --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_printer_fasta_files --input_bed $input_bed --input_fasta $input_fasta --upstream_span $upstream_span --downstream_span $downstream_span --type $type --out $output_dir")
#myjobid_seff_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --open-mode=append --output=$outfile_printer_fasta_files --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_printer_fasta_files >> $outfile_printer_fasta_files")




####################################



eval "$(conda shell.bash hook)"

conda activate my_gimme_scan

reference_TFmotif_collections=$(echo "HOMER,JASPAR2020_vertebrates,CIS-BP")

a=($(echo "$reference_TFmotif_collections" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        reference_TFmotif_collections_sel=${i}
        echo "$reference_TFmotif_collections_sel"

	type=$(echo "$reference_TFmotif_collections_sel""_""gimme_scan""_""$analysis")
	outfile_gimme_scan=$(echo "$Log_files""outfile_""$type"".log")
	touch $outfile_gimme_scan
	echo -n "" > $outfile_gimme_scan
	name_gimme_scan=$(echo "$type""_job")
	seff_name=$(echo "seff""_""$type")

	output_file=$(echo "$output_dir""/""$reference_TFmotif_collections_sel"".bed")
	fasta_file_for_prediction=$(echo "$output_dir""/""ALL_variants_REF_and_ALT.fasta")

	myjobid_gimme_scan=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --job-name=$name_gimme_scan --output=$outfile_gimme_scan --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="gimme scan $fasta_file_for_prediction -p $reference_TFmotif_collections_sel -c 0.85 -n 20 -b > $output_file")
	myjobid_seff_gimme_scan=$(sbatch --dependency=afterany:$myjobid_gimme_scan --open-mode=append --output=$outfile_gimme_scan --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_gimme_scan >> $outfile_gimme_scan")


	if [ $reference_TFmotif_collections_sel = 'CIS-BP' ]; then

	    echo "$reference_TFmotif_collections_sel"" ""needs key to motif assignation"
	    
	    # #### CIS_BP_assignation

     	    type=$(echo "$reference_TFmotif_collections_sel""_""assignation""_""$analysis")
	    outfile_CIS_BP_assignation=$(echo "$Log_files""outfile_""$type"".log")
	    touch $outfile_CIS_BP_assignation
	    echo -n "" > $outfile_CIS_BP_assignation
	    name_CIS_BP_assignation=$(echo "$type""_job")
	    seff_name=$(echo "seff""_""$type")

	    Rscript_CIS_BP_assignation=$(echo "$Rscripts_path""12_CIS_BP_TFmotif_assignation.R")
	    

	    input_file=$(echo "$output_dir""$reference_TFmotif_collections_sel"".bed")
	    TF_key_library=$(echo "/group/soranzo/manuel.tardaguila/Motif_analysis/Reference_files/TF_Information.txt")
	    output_CIS_BP_assigned=$(echo "$output_dir""$reference_TFmotif_collections_sel""_assigned"".bed")
	    TF_Species=$(echo 'Homo_sapiens')
	    TF_Status=$(echo 'D')	    
	    

    	    myjobid_CIS_BP_assignation=$(sbatch --dependency=afterany:$myjobid_gimme_scan --job-name=$name_CIS_BP_assignation --output=$outfile_CIS_BP_assignation --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=256M --parsable --wrap="Rscript $Rscript_CIS_BP_assignation --input_file $input_file --TF_key_library $TF_key_library --output_CIS_BP_assigned $output_CIS_BP_assigned --TF_Species $TF_Species --TF_Status $TF_Status --type $type --out $output_dir")
	    myjobid_seff_CIS_BP_assignation=$(sbatch --dependency=afterany:$myjobid_CIS_BP_assignation --open-mode=append --output=$outfile_CIS_BP_assignation --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_CIS_BP_assignation >> $outfile_CIS_BP_assignation")
	     
	    # myjobid_CIS_BP_assignation=$(sbatch --job-name=$name_CIS_BP_assignation --output=$outfile_CIS_BP_assignation --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=256M --parsable --wrap="Rscript $Rscript_CIS_BP_assignation --input_file $input_file --TF_key_library $TF_key_library --output_CIS_BP_assigned $output_CIS_BP_assigned --TF_Species $TF_Species --TF_Status $TF_Status --type $type --out $output_dir")
	    # myjobid_seff_CIS_BP_assignation=$(sbatch --dependency=afterany:$myjobid_CIS_BP_assignation --open-mode=append --output=$outfile_CIS_BP_assignation --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_CIS_BP_assignation >> $outfile_CIS_BP_assignation")


            echo "->>>$myjobid_CIS_BP_assignation"
	    arr[${#arr[@]}]="$myjobid_CIS_BP_assignation"

	else

	    arr[${#arr[@]}]="$myjobid_gimme_scan"
		
	fi #$reference_TFmotif_collections_sel eq 'CIS-BP'

   

done
 


done_string=$(echo "--dependency=afterany:"""""${arr[@]}"""")
echo "$done_string"

dependency_string=$(echo $done_string|sed -r 's/ /:/g')

echo "$dependency_string"


#### Merge results from the TF motif predictors

eval "$(conda shell.bash hook)"

conda deactivate

type=$(echo "collect_results""_""$analysis")
outfile_collect_TF_motif_results=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_collect_TF_motif_results
echo -n "" > $outfile_collect_TF_motif_results
name_collect_TF_motif_results=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_collect_TF_motif_results=$(echo "$Rscripts_path""28_TFmotif_collection_v2.R")


input_file_HOMER=$(echo "$output_dir""HOMER"".bed")
input_file_JASPAR=$(echo "$output_dir""JASPAR2020_vertebrates"".bed")
input_file_CIS_BP=$(echo "$output_dir""CIS-BP_assigned"".bed")
input_bed=$(echo "$output_dir""VARS.bed")

context_initiator_TFs_string=$(echo "EHF,ELF3,ELF5,ETV4,ETV5,ETS2,ELK1,ETV1,ETV2,ERG,ETS1,ELF1,ELF2,GABPA,HNF1A,HNF1B,IRF3,NF2L2,NF2L1,BATF,PO2F2,PEBB,RUNX2,RUNX1,RUNX3,FOSL1,JUNB,JUND,FOS,FOSB,FOSL2,JUN,BATF3,IRF7,BC11A,SPI1,PU.1,SPIB,STAT2,IRF1,IRF2,IRF9,PRDM1,IRF4,IRF8") ###### https://www.biorxiv.org/content/10.1101/2023.05.05.539543v1.abstract Kribelbauer et al
POST_NOR_subset=$(echo "/group/soranzo/manuel.tardaguila/Bulk_RNA_seq/RITM0022181_2/Time_and_genotype_analysis/POST_NOR_subset.tsv")
DE_edgeR_results=$(echo "/group/soranzo/manuel.tardaguila/Bulk_RNA_seq/RITM0022181_2/Time_and_genotype_analysis/DE_edgeR_results.tsv")

echo "$POST_NOR_subset"

myjobid_collect_TF_motif_results=$(sbatch $dependency_string  --job-name=$name_collect_TF_motif_results --output=$outfile_collect_TF_motif_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_collect_TF_motif_results --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --input_file_CIS_BP $input_file_CIS_BP --input_bed $input_bed --context_initiator_TFs_string $context_initiator_TFs_string --POST_NOR_subset $POST_NOR_subset  --type $type --out $output_dir")
myjobid_seff_collect_TF_motif_results=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --open-mode=append --output=$outfile_collect_TF_motif_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_TF_motif_results >> $outfile_collect_TF_motif_results")

#myjobid_collect_TF_motif_results=$(sbatch --job-name=$name_collect_TF_motif_results --output=$outfile_collect_TF_motif_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_collect_TF_motif_results --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --input_file_CIS_BP $input_file_CIS_BP --input_bed $input_bed --context_initiator_TFs_string $context_initiator_TFs_string --POST_NOR_subset $POST_NOR_subset --type $type --out $output_dir")
#myjobid_seff_collect_TF_motif_results=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --open-mode=append --output=$outfile_collect_TF_motif_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_TF_motif_results >> $outfile_collect_TF_motif_results")



#### TF_motif_stats #############################


type=$(echo "TF_motif_stats""_""$analysis")
outfile_TF_motif_stats=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_TF_motif_stats
echo -n "" > $outfile_TF_motif_stats
name_TF_motif_stats=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_TF_motif_stats=$(echo "$Rscripts_path""29_TFmotif_stats.R")


Table_S6=$(echo "/group/soranzo/manuel.tardaguila/Motif_analysis/""Table_S6.tsv")
TF_motifs_annotated=$(echo "$output_dir""TF_motifs_FINAL_K562_EXP.tsv")
version_string=$(echo 'REF,ALT')
#version_string=$(echo 'REF')
#regulation_string=$(echo 'ATU,Transcriptional_R_ATU,Transcriptional_R')
regulation_string=$(echo 'Transcriptional_R_ATU,Transcriptional_R')
#regulation_string=$(echo 'Transcriptional_R')
Intersect_SNP=$(echo 'YES')
context_initiator_TFs_subset=$(echo "IRF3,IRF7,IRF1,IRF2,IRF9,IRF4,IRF8,RUNX2,RUNX1,RUNX3,FOSL1,FOS,FOSB,FOSL2,PU.1") ###### https://www.biorxiv.org/content/10.1101/2023.05.05.539543v1.abstract Kribelbauer et al

context_only_TFs_subset=$(echo "FOXA1,FOXA2,FOXA3,FOXC,FOXP1,FOXP2,GATA3,MEF2C,LHX1,NKX6-1")  ###### https://www.biorxiv.org/content/10.1101/2023.05.05.539543v1.abstract Kribelbauer et al


 myjobid_TF_motif_stats=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --job-name=$name_TF_motif_stats --output=$outfile_TF_motif_stats --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_TF_motif_stats --Table_S6 $Table_S6 --version_string $version_string --regulation_string $regulation_string --Intersect_SNP $Intersect_SNP --type $type --TF_motifs_annotated $TF_motifs_annotated --context_only_TFs_subset $context_only_TFs_subset --context_initiator_TFs_subset $context_initiator_TFs_subset --out $output_dir")
 myjobid_seff_TF_motif_stats=$(sbatch --dependency=afterany:$myjobid_TF_motif_stats --open-mode=append --output=$outfile_TF_motif_stats --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_TF_motif_stats >> $outfile_TF_motif_stats")

# myjobid_TF_motif_stats=$(sbatch --job-name=$name_TF_motif_stats --output=$outfile_TF_motif_stats --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_TF_motif_stats --Table_S6 $Table_S6 --version_string $version_string --regulation_string $regulation_string --Intersect_SNP $Intersect_SNP --type $type --TF_motifs_annotated $TF_motifs_annotated --context_only_TFs_subset $context_only_TFs_subset --context_initiator_TFs_subset $context_initiator_TFs_subset --out $output_dir")
# myjobid_seff_TF_motif_stats=$(sbatch --dependency=afterany:$myjobid_TF_motif_stats --open-mode=append --output=$outfile_TF_motif_stats --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_TF_motif_stats >> $outfile_TF_motif_stats")



