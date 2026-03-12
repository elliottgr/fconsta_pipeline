
#' Process and filter BLAST hits and generate presence/absence matrices
#'
#' @param input Character. Path to a tab-separated BLAST output file.
#'
#' @param min_pident Numeric. Minimum percent identity required for a BLAST hit
#' to be retained. Hits with `pident` lower than this threshold are discarded
#' before any coverage filtering.
#'
#' @param min_scov Numeric. Minimum subject coverage threshold (percentage).
#' The interpretation of this threshold depends on the value of `allow_gaps`.
#'
#' @param allow_gaps Logical. Controls how subject coverage is calculated and
#' therefore how fragmented or split BLAST hits are handled.
#'
#' If `allow_gaps = TRUE`, the function allows fragmented alignments of the same
#' gene to be combined. In this case, coverage filtering is performed using
#' `total_subject_coverage`, which is calculated by summing the aligned lengths
#' of all BLAST high-scoring segment pairs (HSPs) for the same `(query, gene)`
#' combination. This approach allows genes that are split across multiple
#' alignments (for example due to frameshifts, assembly fragmentation, or
#' insertion/deletion events) to still pass the coverage filter if the combined
#' aligned regions cover at least `min_scov` percent of the subject gene.
#'
#' If `allow_gaps = FALSE`, each alignment fragment is evaluated independently.
#' Coverage filtering is then performed using `simple_subject_coverage`, which
#' represents the coverage of a single BLAST alignment relative to the subject
#' length. In this mode, fragmented hits cannot be combined, meaning that only
#' individual alignments that cover at least `min_scov` percent of the subject
#' gene are retained. This provides a stricter filter and helps avoid detecting
#' genes that are highly fragmented or partially matched.
#'
#' In summary:
#' - `allow_gaps = TRUE`: multiple fragments can be combined to reach the
#'   coverage threshold (more permissive; tolerant to split genes).
#' - `allow_gaps = FALSE`: only single continuous alignments are considered
#'   (more stringent; requires a long contiguous match).
#'
#' @param annotation_file Optional character path to an Excel file containing
#' gene annotations. The file should include columns `Gene`, `Bacteriocin`,
#' and `Predicted functionality`. If provided, these annotations are merged
#' into the BLAST results and the presence/absence matrix.
#'
#' @return A list containing:
#' \describe{
#'   \item{df_ranked}{Filtered and ranked BLAST hits with coverage metrics.}
#'   \item{df_pa_raw}{Presence/absence matrix of genes per genome without annotation.}
#'   \item{df_pa}{Annotated presence/absence matrix (if `annotation_file` provided).}
#'   \item{df_filtered}{Intermediate filtered BLAST hits used for matrix generation.}
#' }
#'
#' @details
#' The function processes BLAST results to identify gene presence across genomes.
#' It calculates two coverage metrics:
#'
#' * `simple_subject_coverage`: coverage from a single BLAST alignment.
#' * `total_subject_coverage`: cumulative coverage across multiple alignments
#' for the same gene and query.
#'
#' Depending on the value of `allow_gaps`, the function uses one of these
#' coverage metrics to determine whether a hit satisfies the coverage threshold.
#'
#' @examples
#' \dontrun{
#' results <- process_blast_hits(
#'   input = "blast_results.tsv",
#'   min_pident = 70,
#'   min_scov = 50,
#'   allow_gaps = TRUE
#' )
#' }

process_blast_hits <- function(
    input = "outputs/blast_outputs.tsv",
    min_pident = 70,
    min_scov = 50,          # Minimum Subject Coverage (e.g., 50%)
    allow_gaps = TRUE,      # TRUE filters by total_coverage, FALSE filters by simple_coverage
    annotation_file = NULL,
    collapse_variants = FALSE # when TRUE, all gene variants get collapsed to the same gene. Otherwise, we generate the PA table while respecting variants
) {
  
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(readxl)
  
  # STEP 1: Read and preprocess BLAST hits with BOTH coverage types
  df_filtered <- read_tsv(input, show_col_types = FALSE) %>%
    dplyr::mutate(
      qid = str_extract(qseqid, "^[^|]+"),
      gene_original = case_when(collapse_variants = FALSE ~ sseqid, collapse_variants = TRUE ~ str_extract(sseqid, "^[^|]+") %>% str_remove_all("[^[:alnum:]]+")),
      #gene_original = str_extract(sseqid, "^[^|]+") %>% str_remove_all("[^[:alnum:]]+"),
      gene_lower = tolower(gene_original),
      #gene_original = sseqid,
      #gene_lower = tolower(sseqid),
      
      # 1. Simple Coverage: Coverage of this single specific alignment (HSP)
      simple_subject_coverage = ((abs(send - sstart) + 1) / slen) * 100
    ) %>%
    filter(pident >= min_pident) %>%    # Filter by percent identity first
    
    # 2. Total Coverage: Combine split genes / frameshifts
    group_by(qid, gene_lower) %>%
    mutate(
      total_subject_aligned = sum(abs(send - sstart) + 1),
      total_subject_coverage = (total_subject_aligned / first(slen)) * 100,
      
      # Cap total coverage at 100% (in case split fragments overlap slightly)
      total_subject_coverage = ifelse(total_subject_coverage > 100, 100, total_subject_coverage)
    ) %>%
    ungroup() %>%
    # 3. Apply the dynamic coverage filter based on allow_gaps parameter
    filter(
      if (allow_gaps) {
        total_subject_coverage >= min_scov
      } else {
        simple_subject_coverage >= min_scov
      }
    )
  
  # STEP 2: Rank best BLAST hits per query-gene pair
  df_ranked <- df_filtered %>%
    dplyr::group_by(qid, gene_lower) %>%
    # Keep the best scoring fragment for the final report
    slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::arrange(desc(bitscore)) %>%
    dplyr::mutate(rank = row_number()) %>%
    # Select both coverage metrics to output so you can compare them!
    dplyr::select(rank, qid, gene_original, gene_lower, bitscore, pident, 
                  simple_subject_coverage, total_subject_coverage, 
                  length, slen, evalue, stitle) %>%
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
  
  # write_tsv(df_ranked, output)
  
  # STEP 4: Create PA matrix using lowercase genes for consistency
  df_pa_raw <- df_filtered %>%
    dplyr::distinct(qid, gene_lower) %>%
    dplyr::mutate(present = 1) %>%
    pivot_wider(names_from = gene_lower, values_from = present, values_fill = 0) %>%
    dplyr::rename(GCA_ID = qid)
  
  # write_tsv(df_pa_raw, pa_output)
  
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
  
  # readr::write_excel_csv(df_pa, pa_annot_output)
  
  return(list(
    df_ranked = df_ranked,
    df_pa_raw = df_pa_raw,
    df_pa = df_pa,
    df_filtered = df_filtered
  ))
}  

