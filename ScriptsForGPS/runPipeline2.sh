#!/bin/bash --login
#dcsrsoft use 20241118

#SBATCH --job-name ELG_Generate_Table
#SBATCH --time 72:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user elliott.greene@unil.ch

module purge
dcsrsoft use 20241118
module load miniforge3/24.11.3-2
module load r-light
conda_init
conda activate /work/FAC/FBM/DBC/slehtine/elliott/bacteriocins/fconsta_pipeline/fcpipe

helpFunction()
{
   echo ""
   echo "Usage: $0 -d parameterd"
   echo -e "\t-d Download genomes"
   exit 1 # Exit script after printing help
}

parameterd=FALSE

## option handling
while getopts "B:" opt
do
   case "$opt" in
      d ) parameterd="$OPTARG" ;; ## when runnning genome download again
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


echo "Running script 01, downloading the genomes"
## downloading genomes 
# if $parameterd; then
#     bash -l scripts/01_download_NCBI_genomes.sh \
#         --NCBI_accession PRJEB3084 \
#         --include genome,gff3,gbff,protein \
#         --output GPS_WGS_data \
#         --include_file_name_in_header \
#         --clean
# fi

# echo "Running script 02, this is the BLAST step"
# ## we already have the reference sequences, just need to run the blast step
# bash scripts/03_BLAST.sh \
#     --blast_db downloaded_data/blast/bacteriocin_BLAST \
#     --blast_db_reference_sequences downloaded_data/UNIPROT_gene_data/protein_sequences/protein_sequences.fasta \
#     --blast_db_type prot \
#     --blast_db_title 'bacteriocin_BLAST_db' \
#     --query_list db/querylist2.tsv \
#     --blast_search blastx \
#     --evalue 1e-3 \
#     --outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore stitle qseq sseq" \
#     --num_threads 6 \
#     --qcov_hsp_perc 0 \
#     --output_file intermediate_blast_outputs_2.tsv

echo "Running script 05, making a PA table"
## generating a PA table
END=88

for i in $(seq 0 $END); do filename=$(printf "inter_blast_outputs_%02d" $i); 
echo Generating table for $filename;
Rscript scripts/05_make_PA_table.R GPS_Intermediate_$i outputs/$filename.tsv;
done
#   
