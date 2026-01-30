#!/bin/bash

# Function to display help message
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
    echo "  --output_file                   (Optional) File to save BLAST results (default: blast_results.tsv)."
    echo
    echo "Example:"
    echo "  $0 --blast_db my_db --blast_db_reference_sequences ref.fasta.gz --blast_db_type prot --blast_db_title 'My Database' --query_list queries.tsv --blast_search blastx --evalue 1e-5 --outfmt \"6 qseqid sseqid pident length qlen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore\" --qcov_hsp_perc 50 --num_threads 4 --output_file results.tsv"
    exit 0
}

# Default values
EVALUE=1e-5
OUTFMT=6
QUERY_LIST="db/query_list.tsv"
BLAST_DB=""
BLAST_DB_REF=""
BLAST_DB_TYPE=""
BLAST_DB_TITLE=""
BLAST_SEARCH=""
QCOV_HSP_PERC="90"
NUM_THREADS=1
OUTPUT_FILE="blast_results.tsv"  # Default output file

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --blast_db) BLAST_DB="$2"; shift 2 ;;
        --blast_db_reference_sequences) BLAST_DB_REF="$2"; shift 2 ;;
        --blast_db_type) BLAST_DB_TYPE="$2"; shift 2 ;;
        --blast_db_title) BLAST_DB_TITLE="$2"; shift 2 ;;
        --query_list) QUERY_LIST="$2"; shift 2 ;;
        --blast_search) BLAST_SEARCH="$2"; shift 2 ;;
        --evalue) EVALUE="$2"; shift 2 ;;
        --outfmt) OUTFMT="$2"; shift 2 ;;
        --qcov_hsp_perc) QCOV_HSP_PERC="$2"; shift 2 ;;
        --num_threads) NUM_THREADS="$2"; shift 2 ;;
        --output_file) OUTPUT_FILE="$2"; shift 2 ;;
        -h|--help) show_help ;;
        *) echo "Unknown argument: $1"; show_help ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$BLAST_DB" || -z "$BLAST_SEARCH" || -z "$QUERY_LIST" ]]; then
    echo "Error: Missing required arguments."
    show_help
fi

# Create BLAST database if it does not exist
if [[ ! -f "${BLAST_DB}.phr" ]]; then
    if [[ -z "$BLAST_DB_REF" || -z "$BLAST_DB_TYPE" || -z "$BLAST_DB_TITLE" ]]; then
        echo "Error: Missing required arguments to build BLAST database."
        show_help
    fi
    echo "Creating BLAST database at $BLAST_DB..."
    
    CLEAN_BLAST_DB_REF=$(mktemp)  # Temporary file for cleaned FASTA headers

    if [[ "$BLAST_DB_REF" == *.gz ]]; then
        gunzip -c "$BLAST_DB_REF" | sed '/^>/ s/ /_/g' > "$CLEAN_BLAST_DB_REF"
    else
        sed '/^>/ s/ /_/g' "$BLAST_DB_REF" > "$CLEAN_BLAST_DB_REF"
    fi

    makeblastdb -in "$CLEAN_BLAST_DB_REF" -dbtype "$BLAST_DB_TYPE" -out "$BLAST_DB" -title "$BLAST_DB_TITLE"
    
    rm "$CLEAN_BLAST_DB_REF"  # Remove temporary file after database creation
else
    echo "Using existing BLAST database at $BLAST_DB..."
fi

# Read query files from TSV file
QUERY_FILES=()
while IFS=$'\t' read -r FILE_PATH; do
    [[ -n "$FILE_PATH" ]] && QUERY_FILES+=("$FILE_PATH")  # Skip empty lines
done < "$QUERY_LIST"

echo $QUERY_FILES

# Clear (or create) the output file
> "$OUTPUT_FILE"

# If a custom outfmt is provided (i.e., not just a number), add the header line to the output file.
if ! [[ "$OUTFMT" =~ ^[0-9]+$ ]]; then
    echo -e "qseqid\tsseqid\tpident\tlength\tqlen\tmismatch\tgapopen\tgaps\tnident\tqstart\tqend\tsstart\tsend\tevalue\tqcovs\tqcovhsp\tbitscore" >> "$OUTPUT_FILE"
fi


# Run BLAST search for each query file after cleaning headers (replace spaces with underscores)
for QUERY in "${QUERY_FILES[@]}"; do
    echo "Cleaning headers in $QUERY..."
    CLEAN_QUERY=$(mktemp)
    sed '/^>/ s/ /_/g' "$QUERY" > "$CLEAN_QUERY"
    
    echo "Running $BLAST_SEARCH on cleaned query file..."
    $BLAST_SEARCH -query "$CLEAN_QUERY" -db "$BLAST_DB" -evalue "$EVALUE" -outfmt "$OUTFMT" -qcov_hsp_perc "$QCOV_HSP_PERC" -num_threads "$NUM_THREADS" >> "$OUTPUT_FILE"
    
    rm "$CLEAN_QUERY"
done

echo "All BLAST searches completed. Results saved in $OUTPUT_FILE."
