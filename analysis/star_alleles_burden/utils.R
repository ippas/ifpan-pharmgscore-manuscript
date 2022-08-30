library(tidyverse)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# # https://github.com/Bioconductor/bioconductor_docker/issues/22
# BiocManager::install("preprocessCore", configure.args="--disable-threading")
library("preprocessCore")


# limits from dbNSFP database
SCORES_LIMITS <- list(
  "cadd_raw" = c("min" = -6.46, "th" = 0, "max" = 18.3),
  "dann" = c("min" = 0, "th" = 0.5, "max" = 1),
  "fathmm_xf" = c("min" = 0, "th" = 0.5, "max" = 1),
  "provean" = c("min" = -14, "th" = 2.5, "max" = 14),
  "mutation_assessor" = c("min" = -5.17, "th" = 1.935, "max" = 6.49)
)

BURDEN_SCORE_NAMES <- c(
  "cadd_raw", "fathmm_xf", "provean", "mutation_assessor"
)

normalize_quantiles <- function(x, target) {
  #' Quantile normalization vector based upon a specified target distribution
  if (all(is.na(x))) {
    NA
  } else {
    y <- normalize.quantiles.use.target(as.matrix(x), target)
    as.vector(y)
  }
}


pass <- function(x) {
  #' Do nothing
  x
}


phred <- function(x) {
  #' Phred values in a vector x, NAs are left as is.

  percent_rank = dense_rank(-x) / length(unique(x[!is.na(x)]))
  -10 * log10(percent_rank)
}


sum_without_na <- function(x) {
  #' Sum values removing missing values. Empty vector is NA.

  if_else(
    all(is.na(x)),
    as.double(NA),
    sum(x, na.rm = TRUE)
  )
}


n_unique_rows <- function(df, by) {
  #' Count distinct rows determined by column 'by'
  df %>%
    distinct(.data[[by]]) %>%
    nrow()
}


convert_pharmvar_scores <- function(path) {
  #' Read score values, set factor levels, rename columns,
  #' invert PROVEAN, and so on
  df <- read_csv(path) %>%
    replace_na(list(`function` = "other")) %>%
    mutate(`function` = factor(
      `function`,
      levels = c(
        "no function", "decreased function", "normal function",
        "increased function", "function not assigned", "uncertain function",
        "unknown function", "other"
      )
    ))
  levels(df$`function`) <- c(
    "no function", "decreased function", "normal function",
    "increased function", "other", "other", "other", "other"
  )

  df_rs <- df %>%
    filter(str_detect(star_allele, "rs")) %>%
    mutate(star_gene = gene, star_alle = star_allele, star_suba = NA)

  df_star <- df %>%
    filter(!str_detect(star_allele, "rs")) %>%
    separate(
      star_allele,
      c("star_gene", "star_alle", "star_suba"),
      sep = "\\*|\\.",
      fill = "right",
      remove = FALSE
    )

  pharmvar <- bind_rows(df_rs, df_star)

  pharmvar <- pharmvar %>%
    select(
      -rsid, -qual, -filters, -info.VI, -evidence_level,
      -dbnsfp.Ensembl_geneid, -dbnsfp.Ensembl_transcriptid, -dbnsfp.genename,
      -dbnsfp.MutationAssessor_score, -dbnsfp.PROVEAN_score
    ) %>%
    rename(gene_name = gene) %>%
    # sub-staralleles
    mutate(
      is_sub = !is.na(star_suba) | str_detect(star_allele, "^rs"),
      .after = star_allele) %>%
    select(-star_gene, -star_alle, -star_suba) %>%
    mutate(
      dbnsfp.provean = -dbnsfp.provean
    ) %>%
    rename_with(.fn = ~ sub("dbnsfp.", "", .x), .cols = starts_with("dbnsfp"))
}


convert_pharmvar_scores2 <- function(path) {
  #' Read score values, set factor levels, rename columns,
  #' invert PROVEAN, and so on
  df <-
    read_tsv(path, col_types = "cc----ccc-d") %>%
    replace_na(list(`function` = "other")) %>%
    mutate(`function` = factor(
      `function`,
      levels = c(
        "no function", "decreased function", "normal function",
        "increased function", "function not assigned", "uncertain function",
        "unknown function", "other"
      )
    ))
  levels(df$`function`) <- c(
    "no function", "decreased function", "normal function",
    "increased function", "other", "other", "other", "other"
  )

  df_rs <- df %>%
    filter(str_detect(star_allele, "rs")) %>%
    mutate(star_gene = gene, star_alle = star_allele, star_suba = NA)

  df_star <- df %>%
    filter(!str_detect(star_allele, "rs")) %>%
    separate(
      star_allele,
      c("star_gene", "star_alle", "star_suba"),
      sep = "\\*|\\.",
      fill = "right",
      remove = FALSE
    )

  pharmvar <- bind_rows(df_rs, df_star)

  pharmvar <- pharmvar %>%
    rename(gene_name = gene) %>%
    # sub-staralleles
    mutate(
      is_sub = !is.na(star_suba) | str_detect(star_allele, "^rs"),
      .after = star_allele) %>%
    select(-star_gene, -star_alle, -star_suba)
}


convert_healthy_sportwgs_scores <- function(path) {
  #' convert healthy for easy join with pharm var
  df <-
    read_csv(path) %>%
    select(
      -rsid, -a_index, -was_split, -within_gene, -GT
    ) %>%
    unite("star_allele", s, gene_name, remove = FALSE) %>%
    mutate(
      `function` = "healthy",
      is_sub = TRUE,
      dbnsfp.provean = -dbnsfp.provean
    ) %>%
    rename_with(.fn = ~ sub("dbnsfp.", "", .x), .cols = starts_with("dbnsfp"))
}


make_score <- function(gene_names, dist_exon_th = Inf, dist_gene_th = Inf) {
  target_distribution <- exp(seq(-5, 5, length.out = 1000)) - exp(-5)
  # plot(seq(0, 1, length.out = length(target_distribution)), target_distribution)

  gene_scores <- list()
  for (g_name in gene_names){
    g_path <- paste0("data/pharmacogenetic-score/gene-scores/", g_name, ".tsv.bgz")
    if (file.exists(g_path)) {
      gene_scores[[g_name]] <-
        read_tsv(g_path, col_types = "ccddddddldd") %>%
        mutate(provean = -provean) %>%
        mutate(gene_name = g_name) %>%

        mutate(
          across(
            all_of(BURDEN_SCORE_NAMES),
            pass,
            .names = "{.col}_oryg"
          )
        ) %>%

        filter(dist_to_exon <= dist_exon_th) %>%
        filter(dist_to_gene <= dist_gene_th) %>%

        mutate(
          across(
            all_of(BURDEN_SCORE_NAMES),
            ~ normalize_quantiles(.x, target_distribution),
          )
        )
    } else {
      stop(paste0("File: ", g_path, " doesn't exist"))
    }
  }

  all_scores <-
    bind_rows(gene_scores) %>%
    relocate(
      gene_name, dist_to_exon, dist_to_gene, gnomad_nfe_AF, in_gnomad,
      .after = alleles
    )
}


compute_burden <- function(df, score_names) {
  #' Group df by patient and gene; sum scores; average these burden components
  df %>%

    group_by(star_allele, `function`, gene_name) %>%
    summarise(
      n_variants = max(across(all_of(score_names), ~ sum(!is.na(.x)))),
      across(all_of(score_names), sum_without_na),
      .groups = "drop"
    ) %>%

    mutate(
      burden = apply(.[score_names], 1, mean, na.rm = TRUE)
    ) %>%
    select(-all_of(score_names))
}
