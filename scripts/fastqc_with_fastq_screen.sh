#!/bin/bash

#SBATCH -p all
# #SBATCH -D .
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 10
#SBATCH -J fastqc_batch_with_fastq_screen
#SBATCH -o %x%j.out
#SBATCH -e %x%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skmcken@emory.edu

fastqscreen=true
fastqc=true
cleanup=true
threads=10

module load FastQC/0.11.4
source ~/tools/miniconda3/bin/activate fastq_screen

#Runs fastq-screen to test fastqs for contamination
## E.g. human/mouse/ecoli in samples that shouldn't be there

# command to get uniq file ID's
filenames=$(ls NBO12250-57608_sRNA/X202SC20021319-Z01-F001/raw_data/*/*.fq.gz | sort | uniq)
numberfiles=$(ls NBO12250-57608_sRNA/X202SC20021319-Z01-F001/raw_data/*/*.fq.gz | sort | uniq |wc -l)
#cd ../
#Not currently true #above currently excludes I1/I2 files. remove *R* to include them

echo -e "Beginning FastQC/MultiQC batch analysis of "$numberfiles" paired-end files at time:\n`date`\n"

if [ "$fastqscreen" = true ]; then
	echo "FastQ Screen `date`"
	while IFS= read -r a; do ##read list of file IDs one at a time
		fastq_screen --aligner bowtie2 --conf ./fastq_screen.conf --subset 100000 --outdir fastq_screen $a
	done <<< "$filenames" ###this is the list from earlier in variable $filenames
	echo "Done with individual files, generating MultiQC report for Fastq_Screen data `date`"
	cd fastq_screen
	multiqc -n contaminant_multiqc_report -c /home/smckenzie/confs/multiqc_config.yaml .
	cd ../
else
	echo -e "FastQ Screen was set to false.\n"
fi


if [ "$fastqc" = true ]; then
	echo "FastQC"
	while IFS= read -r a; do ##read list of file IDs one at a time
		fastqc -t $threads -o ./ $a
	done <<< "$filenames" ###this is the list from earlier in variable $filenames
else
	echo -e "FastQC was set to false.\n"
fi

#conda deactivate
module unload FastQC/0.11.4

#remove intermediate files and put FastQC files into a separate directory


if [ "$cleanup" = true ]; then
	echo "Running cleanup and multiQC operations."
	mkdir fastqc
	mv *_fastqc* fastqc/ #keeps .zip and .html files from FastQC
	##rm *combined* #deletes the intermediate files
	currentdir=$(pwd)
	echo "Done. All FastQC Html files will be in "$currentdir"/fastqc ."
	#end of original QC parsing. Will do MultiQC now

	cd fastqc
	echo "Running MultiQC as: multiqc ."
	multiqc -n overall_multiqc_report .
	#grab the name and read counts from multiqc's summary of fastqc output
#discard useless .0 at end of read counts and make into txt to be imported to excel
#grep "_R1_" multiqc_data/multiqc_general_stats.txt | awk '{print $1,$6}'| cut -d"_" -f1,5|sed 's/\.0//g'> read_counts.txt
	echo "You can delete the .zip and .html files now if desired."
	mkdir fastqc_files
	mv *_fastqc.zip fastqc_files/
	mv *_fastqc.html fastqc_files/
else
	echo -e "Clean-up was set to false.\n"
fi

echo -e "Done at: `date`"
