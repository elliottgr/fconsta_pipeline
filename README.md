# Strep. Genomes Bacteriocin Genes Screening

**Rationale**: The goal of this pipeline is to systematically screen Streptococcus pneumoniae genomes for bacteriocin-related loci, including toxin precursors, immunity proteins, and associated regulatory genes (e.g., the blp and com operons).

References: https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001060
https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2018.02012/full
https://pmc.ncbi.nlm.nih.gov/articles/PMC4517551/

---

## Installation

**Requirements**: Scripts here are written in a combination of bash, python3, and R. 

This pipeline depends on the [NCBI command line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/).

Script 2 has the following python requirements:
* [Biopython](https://biopython.org/wiki/Download)
* [Numpy](https://numpy.org/install/)
* [Requests](https://requests.readthedocs.io/en/latest/user/install/#install)
* [urllib3](https://pypi.org/project/urllib3/)

Script 4 (Option A) requires gawk:
* [gawk](https://www.gnu.org/software/gawk/)

Generating a presence absence table with R requires:
* [tidyverse](https://tidyverse.org/)

## Repository Structure

```text
.
├── db/ ## (example data)
│   ├── bacteriocin_gene_cluster_mapping.xlsx
│   ├── blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa
│   └── blp_and_others_bacteriocin_and_comABCDEFX_gene_collection_prokka_ready.faa.gz
│   └── gene_names.tsv ## names of genes to download from UniProt 
│   └── query_list.tsv ## the fasta files to search the references against
└── scripts/
    ├── 01_download_NCBI_genomes.sh
    ├── 02_download_UNIPROT_genes.py
    ├── 03_BLAST.sh
    ├── 04_filter_BLAST_results.sh
    ├── 04_filter_blast_results.R
    └── format_faa_for_prokka.sh
```

---

## 1. Genome Data Acquisition

**Goal:** Download WGS data from NCBI repository. Default here downloads genomes from [Croucher et al's Massachusetts dataset](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB2632/). Any other NCBI project accession should be usable.

```bash

bash -l scripts/01_download_NCBI_genomes.sh \
    --NCBI_accession PRJEB2632 \
    --include genome,gff3,gbff,protein \
    --output Croucher_WGS_data \
    --include_file_name_in_header \
    --clean

```

You can also take this time to automatically generate a query list by running 

---

## 2. Reference Protein Preparation


**Important:** Database is already in blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa.

**Goal:** Create a database of bacteriocin genes (UniProt + Manual Additions). 

**Instructions:** This step needs a list of tab-seprated gene names to pull reference sequences from UniProt. An example of this is provided in db/gene_names.tsv

### UniProt Download

```bash

python scripts/02_download_UNIPROT_genes.py \
    --input db/gene_names.tsv \
    --output downloaded_data/UNIPROT_gene_data \
    --taxonomy_id 1313 \
    --min_length 10 \
    --max_length 10000 \
    --clean_fasta_head

```

> **Database:** The sample dataset was curated by Florentin from literature-mined genes (e.g., `sccC`), the `comABCDEFX` operon, and manual curation.

> **Final Reference Path:** `db/blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa`

---

## 3. Method A: BLAST Analysis (Genomic vs Protein)

**Rationale:** Detect small or divergent ORFs (like `cibC`) that exist in the genomic DNA but are missed by prediction tools.

### Step 3A: Run BLASTX

```bash

bash scripts/03_BLAST.sh \
    --blast_db downloaded_data/blast/blp_and_others_and_comABCDEFX \
    --blast_db_reference_sequences downloaded_data/UNIPROT_gene_data/protein_sequences/protein_sequences.fasta \
    --blast_db_type prot \
    --blast_db_title 'blp_and_others_bacteriocin_gene_collection' \
    --query_list db/query_list.tsv \
    --blast_search blastx \
    --evalue 1e-3 \
    --outfmt "6 qseqid sseqid pident length qlen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore" \
    --num_threads 6 \
    --qcov_hsp_perc 90 \
    --output_file blast_outputs.tsv

```

### Step 3B: Filter Results (Bash)

```bash

bash scripts/04_filter_BLAST_results.sh \
    -i outputs/blast_outputs.tsv \
    -p 50 \
    -c 50 \
    -o outputs/clean_blast_outputs.tsv

```

### Step 3B (Alternative): Filter & Rank (R)

Using `scripts/04_filter_blast_results.R` to return a presence / absence table of all the BLAST results

```r
Rscript scripts/04_filter_blast_results.R --input outputs/blast_outputs.tsv 
```

<!--
---

## 4. Method B: Bakta Annotation (Alternative)

**Rationale:** Standardized annotation prioritizing the custom bacteriocin database.
**Important:** First step of Bakta annotation pipeline includes ORF calling -> small ORF (e.g., `cibC` will be missed using this appraoch).

### Format Protein Database

```bash
bash -l scripts/format_faa_for_prokka.sh \
    --bakta \
    --input db/blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa

```

### Run Bakta

```bash
bakta GCA_001082985_WGS/data/GCA_001082985.2_6999_1_19_genomic.fna \
    --db /Volumes/Elements/DB/bakta/db \
    --threads 10 \
    --gram + \
    --keep-contig-headers \
    --output bakta_output \
    --proteins db/blp_and_others_bacteriocin_and_comABCDEFX_gene_collection_prokka_ready.faa.gz \
    --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-crispr --skip-plot

```
-->