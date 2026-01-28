#!/usr/bin/env bash
set -euo pipefail

# Default values
INPUT=""
OUTPUT=""
EXPORT_HEADERS_ONLY=false
UPDATE_FROM_FILE=""
MODE=""  # either "prokka" or "bakta"
PROGNAME=$(basename "$0")

usage() {
    cat <<EOF
Usage: $PROGNAME --input <input_file.faa.gz> [OPTIONS]

Options:
  --input <file>             Input gzip multi-fasta file (required)
  --output <file>            Output file (default: input_basename + _prokka_ready.faa.gz)
  --export_headers_only      Instead of reformatting the whole fasta file, export a TSV file with two columns:
                             original header and new header.
  --update_from_file <file>  Provide a two-column TSV file mapping original header to replacement header.
  --prokka                   Use PROKKA header formatting (default if neither --prokka nor --bakta is given).
  --bakta                    Use BAKTA header formatting.
  --help                     Show this help and exit

Description:
  The script expects a multi-fasta file (gzip compressed) in which header lines are either in the simple
  "gene|SeqID" format or in an annotated format such as:

      >sccA|MF990780.1_prot_AXH01218.1_1 [gene=sccA] [locus_tag=streptococcinC_00001] [protein=lactococcin 972 family bacteriocin SccA] [protein_id=AXH01218.1] [location=1..348] [gbkey=CDS]

  If --update_from_file is NOT provided, the script will parse each header and reformat it as follows:

  For PROKKA mode:
      >SeqID ~~~gene~~~~~~

  For BAKTA mode, two cases are handled:

    1. **Simple header (gene|SeqID):**
         >SeqID gene~~~bacteriocin~~~UniProtKB:SeqID

    2. **Annotated header (contains "[gene=" and "[protein="):**
         >protein_id gene~~~protein~~~GENEBANK:protein_id

  If --update_from_file is provided, then for any header found in the mapping file the replacement header is used.
  Headers not found in the update file are processed normally.

  Note: The flags --prokka and --bakta are mutually exclusive.
EOF
}

# Parse command-line arguments
if [[ "$#" -eq 0 ]]; then
    usage
    exit 1
fi

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input)
            INPUT="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --export_headers_only)
            EXPORT_HEADERS_ONLY=true
            shift
            ;;
        --update_from_file)
            UPDATE_FROM_FILE="$2"
            shift 2
            ;;
        --prokka)
            if [[ -n "$MODE" ]]; then
                echo "Error: Cannot specify both --prokka and --bakta" >&2
                exit 1
            fi
            MODE="prokka"
            shift
            ;;
        --bakta)
            if [[ -n "$MODE" ]]; then
                echo "Error: Cannot specify both --prokka and --bakta" >&2
                exit 1
            fi
            MODE="bakta"
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown parameter: $1" >&2
            usage
            exit 1
            ;;
    esac
done

# Default to PROKKA mode if not specified.
if [[ -z "$MODE" ]]; then
    MODE="prokka"
fi

# Check that required input is set.
if [[ -z "$INPUT" ]]; then
    echo "Error: --input is required" >&2
    usage
    exit 1
fi

# Set default output filename if not provided.
if [[ -z "$OUTPUT" ]]; then
    base=$(basename "$INPUT")
    base="${base%.gz}"
    base="${base%.*}"
    OUTPUT="${base}_prokka_ready.faa.gz"
fi

# Set up update mapping arrays.
update_keys=()
update_vals=()

if [[ -n "${UPDATE_FROM_FILE:-}" ]]; then
    if [[ ! -f "$UPDATE_FROM_FILE" ]]; then
        echo "Error: Update mapping file '$UPDATE_FROM_FILE' not found." >&2
        exit 1
    fi
    while IFS=$'\t' read -r orig new_hdr; do
        [[ -z "$orig" ]] && continue
        update_keys+=("$orig")
        update_vals+=("$new_hdr")
    done < "$UPDATE_FROM_FILE"
fi

# Function to process a header line (input without the leading '>').
process_header() {
    local orig="$1"
    local new_header=""
    local found=0

    # First, check if the header is in the update mapping.
    if [[ ${#update_keys[@]} -gt 0 ]]; then
        for i in "${!update_keys[@]}"; do
            if [[ "${update_keys[$i]}" == "$orig" ]]; then
                new_header="${update_vals[$i]}"
                found=1
                break
            fi
        done
    fi

    if [[ $found -eq 0 ]]; then
        # Distinguish annotated from simple headers.
        if echo "$orig" | grep -q "\[gene=" && echo "$orig" | grep -q "\[protein="; then
            # Annotated header case.
            gene=$(echo "$orig" | sed -n 's/.*\[gene=\([^]]*\)\].*/\1/p')
            protein=$(echo "$orig" | sed -n 's/.*\[protein=\([^]]*\)\].*/\1/p')
            protein_id=$(echo "$orig" | sed -n 's/.*\[protein_id=\([^]]*\)\].*/\1/p')
            # For BAKTA mode with annotated header:
            new_header="${protein_id} ${gene}~~~${protein}~~~GENEBANK:${protein_id}"
        else
            # Assume a simple header in the form "gene|SeqID"
            IFS="|" read -r gene seqid <<< "$orig"
            if [[ -z "$seqid" ]]; then
                seqid="$gene"
                gene=""
            fi
            if [[ "$MODE" == "bakta" ]]; then
                new_header="${seqid} ${gene}~~~bacteriocin~~~UniProtKB:${seqid}"
            else
                new_header="${seqid} ~~~${gene}~~~~~~"
            fi
        fi
    fi

    echo "$new_header"
}

# Determine which command to use to uncompress the file.
if command -v gzcat >/dev/null 2>&1; then
    GUNZIP_CMD="gzcat"
elif command -v zcat >/dev/null 2>&1; then
    GUNZIP_CMD="zcat"
else
    echo "Error: Neither gzcat nor zcat is available." >&2
    exit 1
fi

# Create a temporary file for output.
tmp_out=$(mktemp)

# Process the input file.
if $EXPORT_HEADERS_ONLY; then
    $GUNZIP_CMD "$INPUT" | while read -r line; do
        if [[ "$line" =~ ^\> ]]; then
            orig_header="${line#>}"
            new_hdr=$(process_header "$orig_header")
            echo -e "${orig_header}\t${new_hdr}"
        fi
    done > "$tmp_out"
else
    $GUNZIP_CMD "$INPUT" | while read -r line; do
        if [[ "$line" =~ ^\> ]]; then
            orig_header="${line#>}"
            new_hdr=$(process_header "$orig_header")
            echo ">$new_hdr"
        else
            echo "$line"
        fi
    done > "$tmp_out"
fi

# Compress the output into the target file.
gzip -c "$tmp_out" > "$OUTPUT"
rm "$tmp_out"

echo "Output written to: $OUTPUT"
