# Strep. Genomes Bacteriocin Genes Screening

**Rationale**: The goal of this pipeline is to systematically screen Streptococcus pneumoniae genomes for bacteriocin-related loci, including toxin precursors, immunity proteins, and associated regulatory genes (e.g., the blp and com operons).

References: https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001060
https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2018.02012/full
https://pmc.ncbi.nlm.nih.gov/articles/PMC4517551/

---

## Repository Structure

```text
.
├── db/
│   ├── bacteriocin_gene_cluster_mapping.xlsx
│   ├── blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa
│   └── blp_and_others_bacteriocin_and_comABCDEFX_gene_collection_prokka_ready.faa.gz
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

**Goal:** Download WGS WGS data from any NCBI repository.
```bash
bash -l scripts/01_download_NCBI_genomes.sh \
    --NCBI_accession PRJEB2632 \
    --include genome,gff3,gbff,protein \
    --output 01_PRJEB2632_WGS \
    --include_file_name_in_header \
    --clean

```

---

## 2. Reference Protein Preparation

**Important:** Database is already in blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa.

*Goal:* Create a database of bacteriocin genes (UniProt + Manual Additions).

### UniProt Download

```bash
python scripts/02_download_UNIPROT_genes.py \
    --input gene_names.tsv \
    --output UNIPROT_genes_names_1313 \
    --taxonomy_id 1313 \
    --min_length 10 \
    --max_length 10000 \
    --clean_fasta_head

```

> **Database:** The dataset was supplemented with literature-mined genes (e.g., `sccC`), the `comABCDEFX` operon, and manual curation - see `https://unils-my.sharepoint.com/:x:/g/personal/florentin_constancias_unil_ch/ERcNJYUyCwhDsN0oRiaHzWAB9Y45Zx_NKS0r1Zr7ilXP1w?e=NMbqha` sheet `relevant_genes_2`.
> 
> **Final Reference Path:** `db/blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa`

---

## 3. Method A: BLAST Analysis (Genomic vs Protein)

**Rationale:** Detect small or divergent ORFs (like `cibC`) that exist in the genomic DNA but are missed by prediction tools.

### Step A: Run BLASTX

```bash

bash scripts/03_BLAST.sh \
    --blast_db blast/blp_and_others_and_comABCDEFX \
    --blast_db_reference_sequences db/blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa \
    --blast_db_type prot \
    --blast_db_title 'blp_and_others_bacteriocin_gene_collection' \
    --query_list input_genomic_full.fna \
    --blast_search blastx \
    --evalue 1e-3 \
    --outfmt "6 qseqid sseqid pident length qlen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore" \
    --num_threads 6 \
    --output_file 2026/croutcher_gene_aminoacid_blp_and_others_bacteriocin_input_genomic_fullfna_evalu10e3_hsp0.tsv

```

### Step B: Filter Results (Bash)

```bash
bash scripts/04_filter_BLAST_results.sh \
    -i 2026/croutcher_gene_aminoacid_blp_and_others_bacteriocin_input_genomic_fullfna_evalu10e3_hsp0.tsv \
    -p 50 \
    -c 50 \
    -o 2026/croutcher_gene_aminoacid_blp_and_others_bacteriocin_input_genomic_fullfna_evalu10e3_hsp0_filtered.tsv

```

### Step B Alternative: Filter & Rank (R)

Using `scripts/04_filter_blast_results.R`.


```r
process_blast_hits <- function(
    input, 
    min_pident = 70, 
    min_qcovs = 50,
    colnames_vect = c("qseqid", "sseqid", "pident", "length", "qlen", "mismatch", "gapopen", "gaps",
                      "nident", "qstart", "qend", "sstart", "send", "evalue", "qcovs", "qcovhsp", "bitscore"),
    annotation_file = NULL
) {
  library(dplyr)
  library(readr)
  library(stringr)

  # Load and calculate coverage/clean names
  df_filtered <- read_tsv(input, col_names = colnames_vect, show_col_types = FALSE) %>%
    mutate(
      cov = (nident / length) * 100,
      qid = str_extract(qseqid, "^[^|]+"),
      gene_original = str_extract(sseqid, "(?<=|)[^|]+(?=|)") %>% str_remove_all("[^[:alnum:]]+"),
      gene_lower = tolower(gene_original)
    ) %>%
    filter(pident >= min_pident, cov >= min_qcovs)

  # Rank and Rename using explicit dplyr:: calls to prevent masking errors
  df_ranked <- df_filtered %>%
    group_by(qid, gene_lower) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(bitscore)) %>%
    mutate(rank = row_number()) %>%
    dplyr::select(rank, qid, gene_original, gene_lower, bitscore, pident, qcovs, length, evalue) %>%
    dplyr::rename(GCA_ID = qid, Gene_Name = gene_original)

  return(df_ranked)
}

# Execution
results <- process_blast_hits(input = "2026/croutcher_gene_aminoacid_blp_and_others_bacteriocin_input_genomic_fullfna_evalu10e3_hsp0.tsv")

```

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
