#!/bin/bash

#SBATCH -p all
# #SBATCH -D .
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 10
#SBATCH -J trimmomatic
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skmcken@emory.edu


source activate RNAseq

readType=SE # set to PE or SE for paired-end or single-end processing, respectively
threads=10
fileList=$(ls NBO12250-57608_sRNA/X202SC20021319-Z01-F001/raw_data/*/*.fq.gz | sort | uniq)  #fill in folder path to input files. Can contain wildcards.
outputDir=trimmed_reads
adapterFasta=/home/smckenzie/adapters/full_adapters.fa
trimCommand="ILLUMINACLIP:$adapterFasta:2:30:10 MINLEN:18"
cutadapt=true
cutadaptSeq="AGATCGGAAG"


#Collects files and passes them off to trimmomatic for trimming
mkdir -p $outputDir

while read filepath ; do
    filename=$(echo $filepath | rev | cut -d "/" -f 1 | rev)
    echo "trimming file "$filename
    if [ "$cutadapt" = "true" ] ; then
        cutadapt -a $cutadaptSeq -j $threads -o $outputDir/$filename.cutadapt.fastq $filepath
        trimmomatic $readType -phred33 -threads $threads $outputDir/$filename.cutadapt.fastq $outputDir/$filename.trimmed.fastq $trimCommand
    else
        trimmomatic $readType -phred33 -threads $threads $filepath $outputDir/$filename.trimmed.fastq $trimCommand
    fi
done <<< "$fileList"
