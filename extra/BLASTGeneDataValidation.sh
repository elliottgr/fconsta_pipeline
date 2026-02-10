## The method included here downloads multiple reference sequences per gene
## This file takes these reference sequences and BLASTS them against 
## themselves, which allows us to identify whether our reference sequences
## have any strange annotations. E.g., we can rule out any references where
## they are more closely identified with a different gene than their label
## Most of the code is adapted from script 3, which was written by Florentin (very smart fella). :)

show_help() {
    echo "Usage: $0 --blast_db db_path --blast_db_reference_sequences ref_sequences.fasta[.gz] --blast_db_type db_type --blast_db_title db_title --query_list query_list.tsv --blast_search blastx|blastp|blastn --evalue value --outfmt format --qcov_hsp_perc value --num_threads value [--output_file output.tsv]"
    echo
    echo "Arguments:"
    echo "  --blast_db                      Path where the BLAST database should be stored. If not found, it will be created."
    echo "  --blast_db_reference_sequences  FASTA file used to build the database if it does not exist (can be gzipped)."
    echo "  --blast_db_type                 Type of database (-dbtype flag of makeblastdb)."
    echo "  --blast_db_title                Title for the BLAST database (-title flag of makeblastdb)."
    echo "  --query_list                    Path to a TSV file containing a list of FASTA nucleotide or protein query files."
    echo "  --blast_search                  BLAST search program to use (blastx, blastp, or blastn)."
    echo "  --evalue                        E-value threshold for the BLAST search (-evalue flag)."
    echo "  --outfmt                        Output format for BLAST results (-outfmt flag)."
    echo "  --qcov_hsp_perc                 Query coverage per high-scoring segment percentage (-qcov_hsp_perc flag)."
    echo "  --num_threads                   Number of threads to use for BLAST search (-num_threads flag)."
    echo "  --output_dir                    (Optional) Directory to save within (default: /outputs/)."
    echo "  --output_file                   (Optional) File to save BLAST results (default: blast_results.tsv)."
    echo
    echo "Example with respect to data generated during the tutorial:"
    echo "  $0 --blast_db ValidationDB --blast_db_reference_sequences ../downloaded_data/UNIPROT_gene_data/protein_sequences/protein_sequences.fasta --blast_db_type prot --blast_db_title 'ValidationDB' --query_list queries.tsv --blast_search blastp --evalue 1e-5 --outfmt \"6 qseqid sseqid pident length qlen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore\" --qcov_hsp_perc 50 --num_threads 4 --output_file extra/validation_blast_outputs.tsv"
    exit 0
}
EVALUE=1e-5
OUTFMT=6
QUERY_LIST=""
BLAST_DB=""
BLAST_DB_REF=""
BLAST_DB_TYPE=""
BLAST_DB_TITLE=""
BLAST_SEARCH=""
QCOV_HSP_PERC="50"
NUM_THREADS=6
OUTPUT_DIR="outputs"
OUTPUT_FILE="blast_results.tsv"  # Default output file