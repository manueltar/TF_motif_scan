#!/bin/bash

input_fasta=$1
analysis=$2
MASTER_ROUTE=$3

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0



output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

#rm -rf $output_dir
#mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files

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

	myjobid_gimme_scan=$(sbatch --job-name=$name_gimme_scan --output=$outfile_gimme_scan --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="gimme scan $input_fasta -p $reference_TFmotif_collections_sel -c 0.85 -n 20 -b > $output_file")
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

conda deactivate

type=$(echo "collect_results""_""$analysis")
outfile_collect_TF_motif_results=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_collect_TF_motif_results
echo -n "" > $outfile_collect_TF_motif_results
name_collect_TF_motif_results=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_collect_TF_motif_results=$(echo "$Rscripts_path""13_TFmotif_collection.R")


input_file_HOMER=$(echo "$output_dir""HOMER"".bed")
input_file_JASPAR=$(echo "$output_dir""JASPAR2020_vertebrates"".bed")
input_file_CIS_BP=$(echo "$output_dir""CIS-BP_assigned"".bed")

myjobid_collect_TF_motif_results=$(sbatch $dependency_string --job-name=$name_collect_TF_motif_results --output=$outfile_collect_TF_motif_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_collect_TF_motif_results --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --input_file_CIS_BP $input_file_CIS_BP --type $type --out $output_dir")
myjobid_seff_collect_TF_motif_results=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --open-mode=append --output=$outfile_collect_TF_motif_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_TF_motif_results >> $outfile_collect_TF_motif_results")


#myjobid_collect_TF_motif_results=$(sbatch --job-name=$name_collect_TF_motif_results --output=$outfile_collect_TF_motif_results --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_collect_TF_motif_results --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --input_file_CIS_BP $input_file_CIS_BP --type $type --out $output_dir")
#myjobid_seff_collect_TF_motif_results=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --open-mode=append --output=$outfile_collect_TF_motif_results --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_TF_motif_results >> $outfile_collect_TF_motif_results")



#### cross the results with expression data from K-562

type=$(echo "cross_with_K562_EXP""_""$analysis")
outfile_cross_with_K562_EXP=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_cross_with_K562_EXP
echo -n "" > $outfile_cross_with_K562_EXP
name_cross_with_K562_EXP=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_cross_with_K562_EXP=$(echo "$Rscripts_path""14_TFmotif_cross_with_expression_data.R")


TF_collapsed=$(echo "$output_dir""TF_motifs_collapsed.tsv")
POST_NOR_subset=$(echo "/group/soranzo/manuel.tardaguila/Bulk_RNA_seq/RITM0022181_2/Time_and_genotype_analysis/POST_NOR_subset.tsv")
DE_edgeR_results=$(echo "/group/soranzo/manuel.tardaguila/Bulk_RNA_seq/RITM0022181_2/Time_and_genotype_analysis/DE_edgeR_results.tsv")


myjobid_cross_with_K562_EXP=$(sbatch --dependency=afterany:$myjobid_collect_TF_motif_results --job-name=$name_cross_with_K562_EXP --output=$outfile_cross_with_K562_EXP --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_cross_with_K562_EXP --TF_collapsed $TF_collapsed --POST_NOR_subset $POST_NOR_subset --DE_edgeR_results $DE_edgeR_results --type $type --out $output_dir")
 myjobid_seff_cross_with_K562_EXP=$(sbatch --dependency=afterany:$myjobid_cross_with_K562_EXP --open-mode=append --output=$outfile_cross_with_K562_EXP --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_cross_with_K562_EXP >> $outfile_cross_with_K562_EXP")

#myjobid_cross_with_K562_EXP=$(sbatch --job-name=$name_cross_with_K562_EXP --output=$outfile_cross_with_K562_EXP --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_cross_with_K562_EXP --TF_collapsed $TF_collapsed --POST_NOR_subset $POST_NOR_subset --DE_edgeR_results $DE_edgeR_results --type $type --out $output_dir")
# myjobid_seff_cross_with_K562_EXP=$(sbatch --dependency=afterany:$myjobid_cross_with_K562_EXP --open-mode=append --output=$outfile_cross_with_K562_EXP --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_cross_with_K562_EXP >> $outfile_cross_with_K562_EXP")
	    
#### Merge results from the TF motif predictors


type=$(echo "GR_intersect_and_Sushi_graphs""_""$analysis")
outfile_GR_intersect_and_Sushi_graphs=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_GR_intersect_and_Sushi_graphs
echo -n "" > $outfile_GR_intersect_and_Sushi_graphs
name_GR_intersect_and_Sushi_graphs=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_GR_intersect_and_Sushi_graphs=$(echo "$Rscripts_path""15_GR_SNP_intersection_and_Sushi_plots.R")


input_file_HOMER=$(echo "$output_dir""HOMER"".bed")
input_file_JASPAR=$(echo "$output_dir""JASPAR2020_vertebrates"".bed")
input_file_CIS_BP=$(echo "$output_dir""CIS-BP_assigned"".bed")
input_fasta=$(echo "$input_fasta")
TF_collapsed_EXP_crossed=$(echo "$output_dir""TF_motifs_collapsed_crossed_with_K562_EXP.tsv")
TF_collapsed_EXP_and_DE_crossed=$(echo "$output_dir""TF_motifs_collapsed_crossed_with_K562_EXP_DE_in_ctrl.tsv")
ROI=$(echo '7:101499894-101499973__Del_80bp,7:101499919-101499938__Del_16bp')
SNP=$(echo '7:101499930')

echo "$ROI"

myjobid_GR_intersect_and_Sushi_graphs=$(sbatch --dependency=afterany:$myjobid_cross_with_K562_EXP --job-name=$name_GR_intersect_and_Sushi_graphs --output=$outfile_GR_intersect_and_Sushi_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_GR_intersect_and_Sushi_graphs --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --input_fasta $input_fasta --input_file_CIS_BP $input_file_CIS_BP --TF_collapsed_EXP_crossed $TF_collapsed_EXP_crossed --TF_collapsed_EXP_and_DE_crossed $TF_collapsed_EXP_and_DE_crossed --ROI $ROI --SNP $SNP --type $type --out $output_dir")
myjobid_seff_GR_intersect_and_Sushi_graphs=$(sbatch --dependency=afterany:$myjobid_GR_intersect_and_Sushi_graphs --open-mode=append --output=$outfile_GR_intersect_and_Sushi_graphs --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_GR_intersect_and_Sushi_graphs >> $outfile_GR_intersect_and_Sushi_graphs")


#myjobid_GR_intersect_and_Sushi_graphs=$(sbatch --job-name=$name_GR_intersect_and_Sushi_graphs --output=$outfile_GR_intersect_and_Sushi_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_GR_intersect_and_Sushi_graphs --input_file_HOMER $input_file_HOMER --input_file_JASPAR $input_file_JASPAR --input_fasta $input_fasta --input_file_CIS_BP $input_file_CIS_BP --TF_collapsed_EXP_crossed $TF_collapsed_EXP_crossed --TF_collapsed_EXP_and_DE_crossed $TF_collapsed_EXP_and_DE_crossed --ROI $ROI --SNP $SNP --type $type --out $output_dir")
#myjobid_seff_GR_intersect_and_Sushi_graphs=$(sbatch --dependency=afterany:$myjobid_GR_intersect_and_Sushi_graphs --open-mode=append --output=$outfile_GR_intersect_and_Sushi_graphs --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_GR_intersect_and_Sushi_graphs >> $outfile_GR_intersect_and_Sushi_graphs")

