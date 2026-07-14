#!/bin/bash --login
#dcsrsoft use 20241118

#SBATCH --job-name fixing_files
#SBATCH --time 72:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user elliott.greene@unil.ch
END=88
for i in $(seq 1 $END); do filename=$(printf "inter_blast_outputs_%02d" $i); 
echo Fixing $filename;
sed '1s/^/qseqid\tsseqid\tpident\tlength\tqlen\tslen\tmismatch\tgapopen\tgaps\tnident\tqstart\tqend\tsstart\tsend\tevalue\tqcovs\tqcovhsp\tbitscore\tstitle\tqseq\tsseq\n/' outputs/$filename > outputs/$filename.tsv ;
done