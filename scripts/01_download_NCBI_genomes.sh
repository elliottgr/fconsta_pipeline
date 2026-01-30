#!/bin/bash

# Function to display script usage and example
usage() {
    echo "Usage: $0 --NCBI_accession <accession_number> --include <file_types> --output <output_dir> [--clean] [--include_file_name_in_header]"
    echo "Example: $0 --NCBI_accession GCA_000001735.1 --include genome,gff3,gbff --output my_genomes --clean --include_file_name_in_header"
    exit 1
}

# Determine OS for sed compatibility
if [[ "$OSTYPE" == "darwin"* ]]; then
    SED_INPLACE=("sed" "-i" "")
else
    SED_INPLACE=("sed" "-i")
fi

# Default values for optional flags
CLEAN=false
INCLUDE_FILE_NAME_IN_HEADER=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --NCBI_accession)
            NCBI_ACCESSION=$2
            shift 2
            ;;
        --include)
            INCLUDE_TYPES=$2
            shift 2
            ;;
        --output)
            OUTPUT_DIR=$2
            shift 2
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --include_file_name_in_header)
            INCLUDE_FILE_NAME_IN_HEADER=true
            shift
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$NCBI_ACCESSION" || -z "$INCLUDE_TYPES" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: --NCBI_accession, --include, and --output are required."
    usage
fi

# Activate conda environment with error handling
eval "$(conda shell.bash hook)"
if ! conda activate ncbi_datasets; then
    echo "Error: Failed to activate 'ncbi_datasets' conda environment."
    echo "Ensure the environment exists and try: conda install -c conda-forge ncbi-datasets-cli"
    echo
    echo "Alternatively, you can install the datasets and dataformats programs manually at:"
    echo "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/"
    echo
    read -p "Continue without conda environment? [y]es or [n]o" -n -r REPLY
    echo
    if [[ $REPLY =~ [Nn]$ ]]
        then
            exit 1
        else
            condamode=1
    fi
fi

# if [ -z "$condaless"]
# conda deactivate
# fi

# Create output directories
mkdir -p "$OUTPUT_DIR" || exit 1
FINAL_DATA_DIR="${OUTPUT_DIR}/data"
mkdir -p "$FINAL_DATA_DIR" || exit 1

# Build include flags from comma-separated list
INCLUDE_FLAGS=()
for type in $(tr ',' ' ' <<< "$INCLUDE_TYPES"); do
    INCLUDE_FLAGS+=(--include "$type")
done

# Define paths for intermediate files
ZIP_FILE="${OUTPUT_DIR}/${NCBI_ACCESSION}_bioproject_dataset.zip"
EXTRACT_DIR="${OUTPUT_DIR}/${NCBI_ACCESSION}_bioproject_dataset"

# Download dataset with error checking
if ! datasets download genome accession "$NCBI_ACCESSION" \
    --assembly-source genbank \
    --assembly-version latest \
    "${INCLUDE_FLAGS[@]}" \
    --dehydrated \
    --filename "$ZIP_FILE"; then
    echo "Error: Dataset download failed. Check:"
    echo "1. Internet connection"
    echo "2. NCBI API status (https://www.ncbi.nlm.nih.gov)"
    echo "3. Valid accession number"
    # if [ -z "$condaless"]
    #     conda deactivate
    # fi
    exit 1
fi

# Validate ZIP file
if [[ ! -f "$ZIP_FILE" ]]; then
    echo "Error: Download failed - ZIP file not created."
    conda deactivate
    exit 1
fi

# Unzip and rehydrate dataset
unzip -q "$ZIP_FILE" -d "$EXTRACT_DIR" || { echo "Unzip failed"; exit 1; }
datasets rehydrate --directory "$EXTRACT_DIR" || { echo "Rehydration failed"; exit 1; }

# Process assembly directories
DATA_SUBDIR="$EXTRACT_DIR/ncbi_dataset/data"
if [[ ! -d "$DATA_SUBDIR" ]]; then
    echo "Error: Data directory missing in downloaded package"
    conda deactivate
    exit 1
fi

for assembly_dir in "$DATA_SUBDIR"/*; do
    if [[ -d "$assembly_dir" ]]; then
# Inside the loop where assembly_dir is processed
base_dir=$(basename "$assembly_dir")
# Sanitize directory name but allow dots (e.g., GCA_000001735.1)
safe_base_dir=$(tr -cd '[:alnum:]-_.' <<< "$base_dir" | cut -c 1-50)
        
        for file in "$assembly_dir"/*; do
            if [[ -f "$file" ]]; then
                ext="${file##*.}"
                fname=$(basename "$file")
                
                # Copy with standardized naming
                if [[ "$fname" == GCA_* ]]; then
                    dest="$FINAL_DATA_DIR/$fname"
                else
                    dest="$FINAL_DATA_DIR/${safe_base_dir}_genomic.${ext}"
                fi
                
                cp "$file" "$dest" || exit 1
                
                # Modify FASTA headers if requested
                if $INCLUDE_FILE_NAME_IN_HEADER; then
                    case "$ext" in
                        fna|faa|fasta)
                            if ! "${SED_INPLACE[@]}" -E \
                                "s/^>(.*)/>${safe_base_dir}|\1/; \
                                s/ /_/g;" "$dest"; then
                                echo "Error modifying headers in $dest"
                                exit 1
                            fi
                            ;;
                    esac
                fi
            fi
        done
    fi
done

# Cleanup if requested
if $CLEAN; then
    rm -rf "$ZIP_FILE" "$EXTRACT_DIR" && echo "Cleaned intermediate files"
fi

#conda deactivate
echo "Success! Output files in: $FINAL_DATA_DIR"