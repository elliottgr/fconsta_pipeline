#' Process BLAST results into ranked hits and optional annotated presence/absence matrix
#'
#' @param input Path to BLAST output TSV file
#' @param output Path to output ranked hit table (TSV). If NULL, uses input filename with "_ranked.tsv"
#' @param pa_output Path to presence/absence matrix (TSV). If NULL, uses input filename with "_PA.tsv"
#' @param min_pident Minimum percent identity to retain a hit (default 70)
#' @param min_qcovs Minimum query coverage to retain a hit (default 50)
#' @param colnames Vector of column names for input file (default = BLAST outfmt 6 extended)
#' @param annotation_file Optional path to annotation xlsx file using `;` separator (default NULL)
#'
#' @return A list with two data frames: df_ranked and df_pa df_pa_raw
#'
#' @export

process_blast_hits <- function(
    # input = "/Users/fconsta3/Documents/Sonja/Projects/Spneumo/all_scripts/croutcher_gene_aminoacid_blp_and_others_bacteriocin_input_genomic_fullfna_evalu10e3_hsp0.tsv",
    input = "blast_outputs.tsv",
    output = NULL,
    pa_output = NULL,
    pa_annot_output = NULL,
    min_pident = 70,
    min_qcovs = 50,
    colnames = c("qseqid", "sseqid", "pident", "length", "qlen", "mismatch", "gapopen", "gaps",
                 "nident", "qstart", "qend", "sstart", "send", "evalue", "qcovs", "qcovhsp", "bitscore"),
    annotation_file = NULL
) {
  
  if (is.null(output)) output <- sub("\\.tsv$", "_ranked.tsv", input)
  if (is.null(pa_output)) pa_output <- sub("\\.tsv$", "_PA.tsv", input)
  if (is.null(pa_annot_output)) pa_annot_output <- sub("\\.tsv$", "_PA_annot.tsv", input)
  
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(readxl)
  
  # STEP 1: Read and preprocess BLAST hits
  df_filtered <- read_tsv(input) %>%
    dplyr::mutate(
      cov = (nident / length) * 100,
      qid = str_extract(qseqid, "^[^|]+"),
      gene_original = str_extract(sseqid, "^[^|]+") %>% str_remove_all("[^[:alnum:]]+"),
      gene_lower = tolower(gene_original)
    ) %>%
    filter(pident >= min_pident, cov >= min_qcovs)
  
  # STEP 2: Rank best BLAST hits per query-gene pair
  df_ranked <- df_filtered %>%
    dplyr::group_by(qid, gene_lower) %>%
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::arrange(desc(bitscore)) %>%
    dplyr::mutate(rank = row_number()) %>%
    dplyr::select(rank, qid, gene_original, gene_lower, bitscore, pident, qcovs, length, evalue) %>%
    dplyr::rename(GCA_ID = qid, Gene_Name = gene_original)
  
  # STEP 3: Load and join annotation case-insensitively
  if (!is.null(annotation_file)) {
    annotation <- read_xlsx(annotation_file, sheet = 1) %>%
      dplyr::select(Gene, Bacteriocin, `Predicted functionality`) %>%
      dplyr::mutate(gene_lower = tolower(Gene)) %>%
      dplyr::distinct()
    
    df_ranked <- df_ranked %>%
      left_join(annotation, by = "gene_lower") %>%
      relocate(Bacteriocin, `Predicted functionality`, .after = Gene_Name) %>%
      dplyr::mutate(
        Bacteriocin = replace_na(Bacteriocin, "Not_provided"),
        `Predicted functionality` = replace_na(`Predicted functionality`, "Not_provided")
      ) %>%
      dplyr::select(-gene_lower)
  }
  
  write_tsv(df_ranked, output)
  
  # STEP 4: Create PA matrix using lowercase genes for consistency
  df_pa_raw <- df_filtered %>%
    dplyr::distinct(qid, gene_lower) %>%
    dplyr::mutate(present = 1) %>%
    pivot_wider(names_from = gene_lower, values_from = present, values_fill = 0) %>%
    dplyr::rename(GCA_ID = qid)
  
  write_tsv(df_pa_raw, pa_output)
  
  df_pa <- df_pa_raw
  
  # STEP 5: Annotate PA matrix
  if (!is.null(annotation_file)) {
    # Retain only unique Gene and gene_lower mapping
    annotation_unique <- annotation %>%
      distinct(gene_lower, .keep_all = TRUE)
    
    # Map back the gene_lower column names in df_pa to original Gene names
    gene_map <- setNames(annotation_unique$Gene, annotation_unique$gene_lower)
    
    # Columns in PA matrix to rename (excluding GCA_ID)
    gene_cols <- setdiff(colnames(df_pa), "GCA_ID")
    renamed_cols <- ifelse(gene_cols %in% names(gene_map), gene_map[gene_cols], gene_cols)
    
    colnames(df_pa)[-1] <- renamed_cols
    
    # Add any missing genes from annotation
    missing_genes <- setdiff(annotation_unique$Gene, colnames(df_pa)[-1])
    df_pa[missing_genes] <- NA
    
    # Final gene column order
    gene_order <- c(annotation_unique$Gene, setdiff(colnames(df_pa)[-1], annotation_unique$Gene))
    df_pa <- df_pa %>%
      select(GCA_ID, all_of(gene_order))
    
    # Add annotation rows
    annotation_rows <- tibble(Gene = colnames(df_pa)[-1]) %>%
      left_join(annotation_unique, by = c("Gene" = "Gene")) %>%
      mutate(
        Bacteriocin = replace_na(Bacteriocin, "Not_provided"),
        `Predicted functionality` = replace_na(`Predicted functionality`, "Not_provided")
      )
    
    row_bacteriocin <- c("Bacteriocin", annotation_rows$Bacteriocin)
    row_functionality <- c("Functionality", annotation_rows$`Predicted functionality`)
    
    df_pa <- bind_rows(
      as_tibble_row(setNames(row_bacteriocin, colnames(df_pa))),
      as_tibble_row(setNames(row_functionality, colnames(df_pa))),
      df_pa %>% mutate(across(-GCA_ID, as.character))
    )
  }
  
  readr::write_excel_csv(df_pa, pa_annot_output)
  
  return(list(
    df_ranked = df_ranked,
    df_pa_raw = df_pa_raw,
    df_pa = df_pa,
    df_filtered = df_filtered
  ))
}

process_blast_hits(input = "blast_outputs.tsv",
                   min_pident = 70, min_qcovs = 50) -> results

gff <- results$df_filtered %>%
  mutate(
    seqid  = sseqid,              # or qseqid if you want hits on query
    source = "BLAST",
    type   = "match",
    start  = pmin(sstart, send),
    end    = pmax(sstart, send),
    score  = sprintf("%.2f", pident),      # keep 2‑dec pident
    strand = ifelse(send >= sstart, "+", "-"),
    phase  = ".",
    # build a minimal attributes string
    attributes = paste0("ID=", row_number(),
                        ";Name=", sub("\\|.*", "", sseqid),  # gene part
                        ";Target=", qseqid,
                        " ", qstart, " ", qend,
                        ";evalue=", evalue)
  ) %>%
  # keep only the nine GFF columns, in order
  select(seqid, source, type, start, end, score, strand, phase, attributes)

write_lines("##gff-version 3", "blast_hits.gff3")
write.table(gff,
            file   = "~/blast_hits.gff3",
            quote  = FALSE,
            sep    = "\t",
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE)


results$df_filtered %>% 
  write_tsv("~/testdf_filtered.tsv")

process_blast_hits() -> results


results$df_pa
