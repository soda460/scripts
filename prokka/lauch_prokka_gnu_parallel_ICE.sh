#!/bin/bash

# Date: 28/01/2020

# Author: Jean-Simon Brouard

# Description: This script allow to 'write' individual PROKKA commands in a txt file which can be used by GNU parallel

# Note: This script contains hard links, notably to the The Comprehensive Antibiotic Resistance Database v. 3.03

# Usage lauchProkka.sh (i) full_path_MOB_suite_folder (ii) full_path_PROKKA_folder (iii) SPECIES_INFOS

# A working example : ./lauch_prokka_gnu_parallel_ICE.sh /data/ext4/dataDP/data_hybrid_assemblies/MOB_suite
# /data/ext4/dataDP/data_hybrid_assemblies/macSyFinder_Res13_analysis/Sequences/Annotations SPECIES_INFOS

# Comments: Here are listed the major changes relative to the normal version of lauchProkka:
# 1) There is no infos relative to plasmids since we are looking for ICEs
# 2) Only the longest contig is kept!! (Ideally, one will want complete or near complete genome to perform CORE genome analysis

# Another comment: Do not forget that SPECIES_INFOS is a ***tab-delimited*** file with (i) strain, (ii) genus and (iii) species infos

# SPECIES INFOS
# 3G3	Escherichia coli
# 3G4	Escherichia coli
# etc!

# A final comment: This will not work with relative paths, use absolute paths!

sourceDir=$(pwd)

rm -r -f $2/input_files
rm -r -f $2/output

mkdir -p $2/input_files

printf "\n"

# Creation and filling of 3 arrays
declare -a dna_ids
declare -a genus_array
declare -a species_array

# Make newlines the only separator
IFS=$'\n'       

j=0
for i in $(cat $3); do
	dna_ids[j]=$(echo $i | cut -f 1 | tr -d '\n')
	genus_array[j]=$(echo  $i | cut -f 2 | tr -d '\n')
	species_array[j]=$(echo $i | cut -f 3 | tr -d '\n')
	((j++))
done

# Reinitialization of IFS
IFS=$' \t\n'	

dna_id_count=0

# For all sub folders (DNA_ID)
for i in "${dna_ids[@]}"; do

	echo $i

	# Make a folder to place fasta 'edited' input files
	mkdir -p $2/input_files/$i
	cd $1/$i
	echo "### processing $i $dna_id_count ${genus_array[dna_id_count]} ${species_array[dna_id_count]} ###" 2>&1 | tee -a $2/log.txt
	
	# List all chromosomal sequences
	for chromosome in `ls -1 chromosome.fasta | cut -f 1 -d '.'`; do

		# This custom script will change the fasta header of the files produced by MOBsuite
		fastaReHeader.pl $chromosome.fasta $i $2/input_files/$i

		# !!! We retain only the biggest contig i.e. the first one
		head -n 2 $2/input_files/$i/chromosome.fasta > $2/input_files/$i/largest_chr_contig.fasta		

		# Here we write the full PROKKA command
		echo "prokka --prefix $i " | tr -d '\n' >> $sourceDir/prokka_chr_commands.txt
		echo "--force --addgenes --outdir $2/$i " | tr -d '\n' >> $sourceDir/prokka_chr_commands.txt
		echo "--genus ${genus_array[dna_id_count]} --species ${species_array[dna_id_count]} --strain $i " | tr -d '\n' >> $sourceDir/prokka_chr_commands.txt
		echo "--proteins /data/ext4/dataDP/db/CARD_DB/protein_fasta_protein_homolog_model.fasta --usegenus --evalue 1e-09 $2/input_files/$i/largest_chr_contig.fasta" >> $sourceDir/prokka_chr_commands.txt

	done

	cd $1

	printf "### End of process for $i ###\n\n"
	((dna_id_count++))

done

printf "### Look at a file named prokka_chr_commands in the folder where you lauch this script and use it with GNU parallel! \n\n"
