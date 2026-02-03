# Strep. Genomes Bacteriocin Genes Screening

**Rationale**: The goal of this pipeline is to screen Streptococcus pneumoniae genomes for bacteriocin-related loci, including toxin precursors, immunity proteins, and associated regulatory genes (*e.g.*, the blp and com operons).

**References**: 
https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001060
https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2018.02012/full
https://pmc.ncbi.nlm.nih.gov/articles/PMC4517551/

---

## Repository Structure

```text
.
├── 01_PRJEB2632_WGS
│   └── data
│       ├── GCA_001082985.2_6999_1_19_genomic.fna
│       ├── GCA_001083045.2_7622_3_64_genomic.fna
│       └── GCA_001083265.2_7553_4_47_genomic.fna
├── README.md
├── blast_out.tsv
├── db
│   ├── bacteriocin_gene_cluster_mapping.xlsx
│   └── blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa
├── genomic_fna.tsv
└── scripts
    ├── 01_download_NCBI_genomes.sh
    ├── 02_download_UNIPROT_genes.py
    ├── 03_BLAST.sh
    ├── 04_filter_BLAST_results.sh
    ├── 04_filter_blast_results.R
    └── format_faa_for_prokka.sh


```

---

## 1. Genome Data Download

A few genomes are already included in the repository for you to test.

**Goal:** Download WGS WGS data from any NCBI repository.


```bash
conda create -n ncbi_datasets -c conda-forge ncbi-datasets-cli
```

```bash
conda activate ncbi_datasets

bash -l scripts/01_download_NCBI_genomes.sh \
    --NCBI_accession PRJEB2632 \
    --include genome,gff3,gbff,protein \
    --output 01_PRJEB2632_WGS \
    --include_file_name_in_header \
    --clean

```


---

## 3. Method A: BLAST Analysis (Genomic vs Protein)

**Rationale:** Detect small or divergent ORFs (like `cibC`) that exist in the genomic DNA but are missed by prediction tools.


Generate file with genomic fasta files path to process `genomic_fna.tsv` - one is included in the repo. pointing towards the few included .fna files.

```bash
realpath 01_PRJEB2632_WGS/data/*genomic.fna > genomic_fna.tsv
```

Install blast in a dedicated and fresh conda environment.

```bash
conda create -n blast bioconda::blast -y 
```

### Step A: Run BLASTX

```bash
conda activate blast 

bash scripts/03_BLAST.sh \
                --blast_db blast/blp_and_others_and_comABCDEFX \
                --blast_db_reference_sequences db/blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa \
                --blast_db_type prot  \
                --blast_db_title 'blp_and_others_bacteriocin_gene_collection'  \
                --query_list genomic_fna.tsv \
                --blast_search blastx \
                --evalue 1e-3 \
                --outfmt "6 qseqid sseqid pident length qlen mismatch gapopen gaps nident qstart qend sstart send evalue qcovs qcovhsp bitscore" \
                --qcov_hsp_perc 0 --num_threads 6 \
                --output_file blast_out.tsv
```

We reported all hits with a minimum evalue of `1e-3`. These raw data needs to be filtered. We do can do that in R.


### Step B: Filter & Rank (R)

Using `scripts/04_filter_blast_results.R`.


```r
setwd("~/Documents/GitHub/Streppneumo_bacteriocins/")
source("scripts/04_filter_blast_results.R")

process_blast_hits(input = "blast_out.tsv",
+                    min_pident = 70, min_qcovs = 50,
+                    annotation_file = NULL) -> results$df_pa

> > results$df_pa
> # A tibble: 616 × 130
>    GCA_ID   scab  scac  scaa  cibc  ciba  cibb  slfa  cbpd comfa  comd  come comfc  comc  blph  blpb  blpg
>    <chr>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
>  1 GCA_00…     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
>  2 GCA_00…     0     0     0     1     1     1     1     1     1     1     1     1     1     1     1     1
>  3 GCA_00…     0     0     0     1     1     1     1     1     1     1     1     1     1     1     1     1
>  4 GCA_00…     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1

```

---

### Step B - Legacy Filter Results (Bash)

```bash
bash scripts/04_filter_BLAST_results.sh \
    -i blast_out.tsv \
    -p 50 \
    -c 50 \
    -o blast_out_filtered.tsv

```


---

## S1. Reference Protein Preparation

**Important:** The database is already in blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa.

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
---

> **Database:** The dataset was supplemented with literature-mined genes (e.g., `sccC`), the `comABCDEFX` operon, and manual curation - see `https://unils-my.sharepoint.com/:x:/g/personal/florentin_constancias_unil_ch/ERcNJYUyCwhDsN0oRiaHzWAB9Y45Zx_NKS0r1Zr7ilXP1w?e=NMbqha` sheet `relevant_genes_2`.
> 
> **Final Reference Path:** `db/blp_and_others_bacteriocin_and_comABCDEFX_gene_collection.faa`


## S2. Method B: Bakta Annotation (Alternative)

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
---
