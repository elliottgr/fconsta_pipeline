import argparse
import requests
from Bio import Entrez, SeqIO
import os
import time
import logging
import gzip
import numpy as np
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Set up logging
def setup_logging(log_level):
    logging.basicConfig(level=getattr(logging, log_level),
                        format='%(asctime)s - %(levelname)s - %(message)s')

# Function to calculate average length of sequences for a gene
def calculate_average_length(sequences):
    lengths = [len(seq) for seq in sequences]
    return np.mean(lengths) if lengths else 0

# Function to remove outliers based on sequence length
def remove_length_outliers(sequences, min_length=None, max_length=None, filter_outliers=False):
    if filter_outliers:
        avg_length = calculate_average_length(sequences)
        # Define outlier thresholds (for example, 2 standard deviations away from the mean)
        std_deviation = np.std([len(seq) for seq in sequences])
        lower_threshold = avg_length - 2 * std_deviation
        upper_threshold = avg_length + 2 * std_deviation

        filtered_sequences = [seq for seq in sequences if lower_threshold <= len(seq) <= upper_threshold]
        logging.info(f"Filtered out {len(sequences) - len(filtered_sequences)} sequences as outliers based on length.")
        sequences = filtered_sequences
    return sequences

# Function to download protein sequences from UniProt using the streaming API
def download_protein_sequence(gene_name, taxonomy_id, output_dir, min_length=None, max_length=None,
                              exclude_fields=None, filter_outliers=False, clean_fasta_head=False):
    query = f"(gene_exact:{gene_name}+AND+(taxonomy_id:{taxonomy_id}))"
    
    # Add length filters if provided
    if min_length is not None or max_length is not None:
        if min_length is not None and max_length is not None:
            length_filter = f"+AND+(length:[{min_length}+TO+{max_length}])"
        elif min_length is not None:
            length_filter = f"+AND+(length:[{min_length}+TO+inf])"
        elif max_length is not None:
            length_filter = f"+AND+(length:[0+TO+{max_length}])"
        query += length_filter
    
    # Add exclusion terms to the query
    if exclude_fields:
        for field, value in exclude_fields.items():
            query += f"+NOT+({field}:{value})"
    
    # Construct the URL for UniProt
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query={query}"

    try:
        response = requests.get(url)
        if response.status_code == 200 and response.text:
            sequences = []
            for line in response.text.split("\n"):
                if line.startswith(">"):
                    if clean_fasta_head:
                        # Modify header: gene_name|UniProt_ID (extracted from the second field after splitting on '|')
                        parts = line.split("|")
                        if len(parts) >= 3:
                            uniprot_id = parts[1]
                            header = f">{gene_name}|{uniprot_id}"
                            sequences.append(header)
                        else:
                            sequences.append(line)
                    else:
                        # Keep the original header
                        sequences.append(line)
                else:
                    sequences.append(line)

            # Remove outliers based on sequence length if needed
            sequences = remove_length_outliers(sequences, min_length, max_length, filter_outliers)
            
            # Save the sequences to a FASTA file
            output_file = os.path.join(output_dir, f"{gene_name}_protein.fasta")
            with open(output_file, "w") as f:
                f.write("\n".join(sequences))
            logging.info(f"Downloaded and filtered protein sequence for {gene_name}")
        else:
            logging.warning(f"No protein sequence found for {gene_name} or failed to retrieve from UniProt.")
    except requests.exceptions.RequestException as e:
        logging.error(f"Error downloading protein sequence for {gene_name}: {e}")

# Function to download nucleic acid sequences from NCBI
def download_nucleotide_sequence(gene_name, taxonomy_id, output_dir, ncbi_email, min_length=None, max_length=None,
                                 exclude_fields=None, filter_outliers=False):
    query = f"(gene_exact:{gene_name}+AND+(taxonomy_id:{taxonomy_id}))"
    
    # Add length filters if provided
    if min_length is not None or max_length is not None:
        if min_length is not None and max_length is not None:
            length_filter = f"+AND+(length:[{min_length}+TO+{max_length}])"
        elif min_length is not None:
            length_filter = f"+AND+(length:[{min_length}+TO+inf])"
        elif max_length is not None:
            length_filter = f"+AND+(length:[0+TO+{max_length}])"
        query += length_filter
    
    # Add exclusion terms to the query
    if exclude_fields:
        for field, value in exclude_fields.items():
            query += f"+NOT+({field}:{value})"
    
    # Construct the URL for NCBI Entrez
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term={query}&retmode=xml&email={ncbi_email}"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            result = response.text
            sequences = [seq for seq in result.split("\n") if seq]  # Capture the sequence content lines
            sequences = remove_length_outliers(sequences, min_length, max_length, filter_outliers)
            if sequences:
                output_file = os.path.join(output_dir, f"{gene_name}_nucleotide.fasta")
                with open(output_file, "w") as f:
                    f.write("\n".join(sequences))
                logging.info(f"Downloaded and filtered nucleotide sequence for {gene_name}")
            else:
                logging.warning(f"No nucleotide sequence found for {gene_name} from NCBI.")
        else:
            logging.warning(f"Failed to retrieve nucleotide sequence for {gene_name} from NCBI.")
    except requests.exceptions.RequestException as e:
        logging.error(f"Error downloading nucleotide sequence for {gene_name}: {e}")

# Function to concatenate FASTA files and compress them into a .gz file
def concatenate_and_compress_fasta(output_dir, output_filename):
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith(".fasta")]
    if fasta_files:
        gz_file = os.path.join(output_dir, f"{output_filename}.fasta.gz")
        with gzip.open(gz_file, "wt") as outfile:
            for fasta_file in fasta_files:
                file_path = os.path.join(output_dir, fasta_file)
                with open(file_path, "r") as infile:
                    outfile.write(infile.read())
        logging.info(f"Concatenated and compressed {len(fasta_files)} FASTA files into {output_filename}.fasta.gz")
    else:
        logging.warning(f"No FASTA files found in {output_dir} to concatenate.")

# Main function
def main():
    parser = argparse.ArgumentParser(description="Download protein and/or nucleic acid sequences for a list of genes.")
    parser.add_argument("--input", required=True, help="Path to the input file containing gene names (one per line).")
    parser.add_argument("--output", default=".", help="Path to the output directory where sequences will be saved (default: current directory).")
    parser.add_argument("--taxonomy_id", required=True, help="Taxonomy ID of the organism (e.g., '1313' for Streptococcus pneumoniae).")
    parser.add_argument("--ncbi_email", help="Email for NCBI Entrez (required for nucleic acid sequences).")
    parser.add_argument("--sequence_type", choices=["amino_acid", "nucleic_acid", "both"],
                        default="amino_acid", help="Type of sequences to download: 'amino_acid', 'nucleic_acid', or 'both' (default: 'amino_acid').")
    parser.add_argument("--min_length", type=int, help="Minimum length of the sequences to download.")
    parser.add_argument("--max_length", type=int, help="Maximum length of the sequences to download.")
    parser.add_argument("--log_level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        default="INFO", help="Set the logging level (default: INFO)")
    parser.add_argument("--exclude", help="Comma-separated list of gene names to exclude from processing.")
    parser.add_argument("--exclude_fields", help="Comma-separated list of field:value pairs for exclusion criteria (e.g., 'keyword:Truncated,keyword:Fragment').")
    parser.add_argument("--filter_outliers", action="store_true", help="Whether to filter out outlier sequences based on length.")
    parser.add_argument("--clean_fasta_head", action="store_true",
                        help="If provided, modify protein FASTA headers to use the format 'gene_name|UniProt_ID'.")
    
    args = parser.parse_args()

    if args.sequence_type in ["nucleic_acid", "both"] and not args.ncbi_email:
        raise ValueError("--ncbi_email is required for downloading nucleic acid sequences.")
    
    if args.min_length and args.max_length and args.min_length > args.max_length:
        raise ValueError("min_length cannot be greater than max_length.")

    setup_logging(args.log_level)

    protein_output_dir = os.path.join(args.output, "protein_sequences")
    nucleotide_output_dir = os.path.join(args.output, "nucleotide_sequences")
    if args.sequence_type in ["amino_acid", "both"]:
        os.makedirs(protein_output_dir, exist_ok=True)
    if args.sequence_type in ["nucleic_acid", "both"]:
        os.makedirs(nucleotide_output_dir, exist_ok=True)

    # Set up a requests session with a retry strategy
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    session.mount('http://', HTTPAdapter(max_retries=retries))
    session.mount('https://', HTTPAdapter(max_retries=retries))

    # Read gene names from the input file
    with open(args.input, "r", errors = "replace") as f:
        gene_names = [line.strip() for line in f.readlines() if line.strip()]

    # Remove duplicates and log the number of unique genes
    unique_genes = list(set(gene_names))
    logging.info(f"Provided gene names: {len(gene_names)}, number of unique gene names: {len(unique_genes)}")

    if args.exclude:
        exclude_list = [gene.strip() for gene in args.exclude.split(",")]
        unique_genes = [gene for gene in unique_genes if gene not in exclude_list]
        logging.info(f"Excluded genes: {exclude_list}, remaining genes to process: {len(unique_genes)}")

    if args.exclude_fields:
        exclude_fields_dict = dict(field_value.split(":") for field_value in args.exclude_fields.split(","))
    else:
        exclude_fields_dict = {}

    # Process each gene
    for gene_name in unique_genes:
        if args.sequence_type in ["amino_acid", "both"]:
            download_protein_sequence(
                gene_name, args.taxonomy_id, protein_output_dir,
                args.min_length, args.max_length, exclude_fields_dict,
                args.filter_outliers, clean_fasta_head=args.clean_fasta_head
            )
        if args.sequence_type in ["nucleic_acid", "both"]:
            download_nucleotide_sequence(
                gene_name, args.taxonomy_id, nucleotide_output_dir, args.ncbi_email,
                args.min_length, args.max_length, exclude_fields_dict, args.filter_outliers
            )

    # Concatenate and compress FASTA files for protein and nucleotide sequences
    if args.sequence_type in ["amino_acid", "both"]:
        concatenate_and_compress_fasta(protein_output_dir, "protein_sequences")
    if args.sequence_type in ["nucleic_acid", "both"]:
        concatenate_and_compress_fasta(nucleotide_output_dir, "nucleotide_sequences")

if __name__ == "__main__":
    main()
