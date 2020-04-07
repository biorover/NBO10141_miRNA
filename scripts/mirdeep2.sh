#!/bin/bash

#SBATCH -p all
# #SBATCH -D .
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 2
#SBATCH -J slurmMiRDeep2
#SBATCH -o %x%j.out
#SBATCH -e %x%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skmcken@emory.edu

#
# Script for analyzing trimmed and QC-ed small RNAseq reads for microRNA identification with miRDeep2 and microRNA quantification via bowtie2 alignment
#
#


#####
#
# Setting variables used later in the script
#
#####

source ~/tools/mirdeep2/bin/activate 	# path to file setting environmental variables for mirdeep2, 
                                     	# e.g. putting mirdeep2/bin into PATH variable and setting perl library locations


bowtieIndex=$(readlink -f /home/smckenzie/genomes/pig/bowtie1/pigbowtie)
genomeFasta=$(readlink -f /home/smckenzie/genomes/pig/genome_cleanNames.fa)
threads=10
readFiles=$(ls trimmed_reads/*.trimmed.fastq)
tempFileName=$(echo /tmp/mirdeep2pipeline$(date) | tr -d " " | tr d ":")
outDir="mirdeep_out"
baseDir=$(pwd)
species="Sus_scrofa" # requires some identifier that will be in all conspecific reference fasta entries

extraMapperOptions="-e -j -m -h" 	# set flag for input format (-e for fastq, -c for fasta), set -m to collapse reads,
                              		# set -j to autopurge non-canonical bases (non ACTGUN), set -k "ADAPTER_SEQ" to clip adapters,
                              		# set -l INT to remove seqs smaller than INT, and set -h to parse illumina fastq to fasta (required if input is fastq)

# For miRDeep2
matureRef=$(readlink -f databaseFiles/Sus_scrofa_mature_mirbase.fa)
matureRefOutgroup=$(readlink -f databaseFiles/Homo_sapiens_mature_mirbase.fa)
precursorsRef=$(readlink -f databaseFiles/Sus_scrofa_hairpin_mirbase.fa)

#####
#
# Running mapper to map reads to genome for ID of known and novel miRNAs
#
#####
if [ 0 -eq 1 ] ; then #can be used to switch off if already run

mkdir -p $outDir
cat $readFiles > $tempFileName
command="mapper.pl $tempFileName $extraMapperOptions -p $bowtieIndex -s $outDir/reads_collapsed.fa -t $outDir/reads_collapsed_vs_genome.arf -v"
echo "running $command"
$command

rm $tempFileName

fi
#####
#
# Running miRDeep2 to identify known and novel miRNAs
#
#####
if [ 0 -eq 1 ] ; then #can be used to switch off if already run

cd $outDir
command="miRDeep2.pl reads_collapsed.fa $genomeFasta reads_collapsed_vs_genome.arf $matureRef $matureRefOutgroup $precursorsRef 2> mirDeep2report.log"
echo "running $command"
$command
cd ../

fi
#####
#
# Mapping reads against miRDeep2 predicted microRNAs
#
#####

source activate miRNAanalysis
cd $outDir
mkdir -p MIR_quantification
cd MIR_quantification
#cat ../mirna_results*/*pres*.fa > precursors.fa
#bowtie2-build precursors.fa precursors
#while read fastq ; do
#    fileroot=$(echo $fastq | rev | cut -d "/" -f 1 | rev | sed -e 's/\.fastq$//' -e 's/\.fq$//' )
#    bowtie2 -x precursors -U $baseDir/$fastq -p $threads | samtools sort -@ $threads > $fileroot.bam
#    samtools view -q 30 -F 2304 $fileroot.bam | awk -v sample=$fileroot '{a[$3] += 1}END{for (i in a) print sample "\t" i "\t" a[i]}' | sort -k2,2n > $fileroot.sample_counts
#done <<< "$readFiles"

cat *.sample_counts | awk '{OFS="\t"}
NR>1{
    map[$2,$1] = $3
    name[$2]++
    value[$1]++
}
END{
    printf "Gene"
    n = asorti(value, v_s)
    for(i=1; i<=n; i++) {
        printf "%s%s", FS, v_s[i]
    }
    print ""
    m = asorti(name, n_s)
    for(i=1; i<=m; i++) { 
        printf "%s", n_s[i]
        for(j=1; j<=n; j++) { 
            printf "%s%s", FS, map[n_s[i],v_s[j]]
        }
        print ""
    }
}' | tr " " "\t" | sed -e 's/\t\t/\t0\t/g' -e 's/\t\t/\t0\t/g' -e 's/\t$/\t0\t/g'  > all_samples.counts
echo "sed " | tr -d "\n" > sedcom.sh
grep "$species" ../result*.csv | awk '{print "-e ;s/"$1"\\t/"$12"\\t/; "}' | tr ";" "'" | tr -d "\n" >> sedcom.sh
echo "all_samples.counts" >> sedcom.sh
bash sedcom.sh > all_samples.named.counts
