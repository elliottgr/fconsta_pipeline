#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo
    echo "Filters BLAST results and reports top 3 hits with consensus information"
    echo
    echo "Options:"
    echo "  -i, --input FILE          Input TSV file (required)"
    echo "  -o, --output FILE         Output file (default: <input>_filtered.tsv)"
    echo "  -p, --min_pident FLOAT    Minimum percent identity (default: 0)"
    echo "  -c, --min_qcovs FLOAT     Minimum query coverage (default: 0)"
    echo "  -n, --colnames STRING     Column names (default: \"qseqid sseqid pident length qlen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore\")"
    echo "  -h, --help                Show this help message"
    echo
    echo "Example:"
    echo "  $0 -i blast_results.tsv -p 30 -c 80 -o filtered_hits.tsv"
    exit 0
}

# Initialize variables with defaults
INPUT_FILE=""
OUTPUT_FILE=""
MIN_PIDENT=0
MIN_QCOVS=0
COL_NAMES="qseqid sseqid pident length qlen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore"

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_FILE="$2"; shift 2 ;;
        -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
        -p|--min_pident) MIN_PIDENT="$2"; shift 2 ;;
        -c|--min_qcovs) MIN_QCOVS="$2"; shift 2 ;;
        -n|--colnames) COL_NAMES="$2"; shift 2 ;;
        -h|--help) show_help ;;
        *) echo "Unknown parameter: $1"; show_help ;;
    esac
done

# Validate required input
if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: Input file is required"
    show_help
fi

# Set default output filename if not provided
if [[ -z "$OUTPUT_FILE" ]]; then
    OUTPUT_FILE="${INPUT_FILE%.*}_filtered.tsv"
fi

# Create array of column names
IFS=' ' read -ra COL_ARRAY <<< "$COL_NAMES"

# Function to get column index by name
get_col_index() {
    local col_name="$1"
    for i in "${!COL_ARRAY[@]}"; do
        if [[ "${COL_ARRAY[$i]}" == "$col_name" ]]; then
            echo $((i + 1))  # AWK uses 1-based indexing
            return
        fi
    done
    echo "Error: Column '$col_name' not found in column names" >&2
    exit 1
}

# Get required column indices
QSEQID_COL=$(get_col_index "qseqid")
SSEQID_COL=$(get_col_index "sseqid")
PIDENT_COL=$(get_col_index "pident")
QCOVS_COL=$(get_col_index "qcovs")

# Process the file using gawk (GNU awk) with composite keys
gawk -F '\t' -v MIN_PIDENT="$MIN_PIDENT" -v MIN_QCOVS="$MIN_QCOVS" \
    -v QSEQID_COL="$QSEQID_COL" -v SSEQID_COL="$SSEQID_COL" \
    -v PIDENT_COL="$PIDENT_COL" -v QCOVS_COL="$QCOVS_COL" '
BEGIN {
    OFS = "\t";
    # Print new header with five columns
    print "GCA_ID", "Full_Header", "Status", "Gene_Name", "Details";
}

NR == 1 { next; }  # Skip header

# Filter rows by thresholds
$PIDENT_COL >= MIN_PIDENT && $QCOVS_COL >= MIN_QCOVS {
    query = $QSEQID_COL;
    
    # Store the full sseqid header for use later
    full_sseqid = $SSEQID_COL;
    
    # Parse sseqid to extract gene name and accession.
    # Here we assume that sseqid is formatted as "gene|accession" (if not, adjust accordingly)
    split($SSEQID_COL, parts, "\\|");
    gene = parts[1];
    accession = (parts[2] ? parts[2] : "unknown");
    
    pident = $PIDENT_COL;
    qcovs = $QCOVS_COL;
    
    # Record query order if first seen
    if (!(query in query_count)) {
        query_order[++total_queries] = query;
    }
    
    # Increase hit count and store data using composite key
    idx = ++query_count[query];
    # For details, build a string with gene, accession, percent identity and query coverage
    hit_str = gene "|" accession "_id" sprintf("%.2f", pident) "_cov" qcovs;
    hits[query, idx] = hit_str;
    genes[query, idx] = gene;
    pidents[query, idx] = pident;
    qcovs_values[query, idx] = qcovs;
    # Save the full sseqid header (useful for constructing Full_Header later)
    headers[query, idx] = full_sseqid;
}

END {
    # Process each query in order
    for (q = 1; q <= total_queries; q++) {
        query = query_order[q];
        hit_count = query_count[query];
        
        # Build sort keys for hits: higher percent identity and query coverage sort first.
        # Here, we create keys that sort in ascending order so that lower keys are better.
        for (i = 1; i <= hit_count; i++) {
            sortkeys[i] = sprintf("%06.2f%03d", 100 - pidents[query, i], 100 - qcovs_values[query, i]);
        }
        
        # Selection sort: sort only top three hits.
        max_hits = (hit_count < 3 ? hit_count : 3);
        for (i = 1; i <= max_hits; i++) {
            min_idx = i;
            for (j = i + 1; j <= hit_count; j++) {
                if (sortkeys[j] < sortkeys[min_idx]) {
                    min_idx = j;
                }
            }
            if (min_idx != i) {
                temp = hits[query, i];
                hits[query, i] = hits[query, min_idx];
                hits[query, min_idx] = temp;
                
                temp = genes[query, i];
                genes[query, i] = genes[query, min_idx];
                genes[query, min_idx] = temp;
                
                temp = headers[query, i];
                headers[query, i] = headers[query, min_idx];
                headers[query, min_idx] = temp;
                
                temp = sortkeys[i];
                sortkeys[i] = sortkeys[min_idx];
                sortkeys[min_idx] = temp;
            }
        }
        
        # Build Details string from the top hits
        details = "";
        delete consensus_check;
        for (i = 1; i <= max_hits; i++) {
            hit = hits[query, i];
            gene_name = genes[query, i];
            details = (details ? details " ; " : "") hit;
            consensus_check[gene_name]++;
        }
        
        # Determine consensus status based on the diversity of gene names
        unique_genes = 0;
        for (k in consensus_check) {
            unique_genes++;
        }
        if (unique_genes == 1) {
            status = "Consensus";
        } else if (unique_genes < max_hits) {
            status = "Partial";
        } else {
            status = "Divergent";
        }
        
        # Construct Full_Header using the query and the top hit sseqid.
        # For example: "GCA_001082985.2|CTM98421.1_putative_membrane_protein_[Streptococcus_pneumoniae]"
        top_full_sseqid = headers[query, 1];
        full_header = query "|" top_full_sseqid;
        
        # Use the gene name from the top hit as Gene_Name
        top_gene = genes[query, 1];
        
        print query, full_header, status, top_gene, details;
    }
}' "$INPUT_FILE" | sort -t $'\t' -k3,3 -k1,1 > "$OUTPUT_FILE"

echo "Filtered results saved to: $OUTPUT_FILE"
