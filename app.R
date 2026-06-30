####################### Sylph MPA Metagenomics Explorer #######################
# GitHub-ready version with package checks, safer validation, reproducibility
# outputs, citation helper, and publication-friendly exports.

# -----------------------------
# 0. Package setup
# -----------------------------
# Set AUTO_INSTALL_PACKAGES <- FALSE if you prefer to install packages manually.
AUTO_INSTALL_PACKAGES <- TRUE

cran_packages <- c(
  "shiny", "bslib", "tidyverse", "vegan", "DT", "plotly", "pals",
  "rmarkdown", "dendextend", "forcats", "gridExtra", "MASS", "htmlwidgets"
)
bioc_packages <- c("phyloseq")

install_missing_packages <- function(cran_pkgs, bioc_pkgs, auto_install = TRUE) {
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (!auto_install && (length(missing_cran) > 0 || length(missing_bioc) > 0)) {
    stop(
      "Missing packages: ",
      paste(c(missing_cran, missing_bioc), collapse = ", "),
      "\nInstall them manually or set AUTO_INSTALL_PACKAGES <- TRUE."
    )
  }

  if (length(missing_cran) > 0) {
    message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }

  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

install_missing_packages(cran_packages, bioc_packages, AUTO_INSTALL_PACKAGES)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(tidyverse)
  library(phyloseq)
  library(vegan)
  library(DT)
  library(plotly)
  library(pals)
  library(rmarkdown)
  library(dendextend)
  library(forcats)
  library(gridExtra)
  library(MASS)
  library(htmlwidgets)
})

options(shiny.maxRequestSize = 500 * 1024^2)
options(width = 120)
APP_VERSION <- "0.3.0-github-ready"
APP_NAME <- "Sylph MPA Metagenomics Explorer"

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  titlePanel("Sylph MPA Metagenomics Explorer"),

  sidebarLayout(
    sidebarPanel(
      h4("1. Data Input"),
      fileInput("files", "Upload Advanced Sylph MPA Files", multiple = TRUE, accept = c(".sylphmpa", ".txt", ".tsv")),
      fileInput("metadata_file", "Optional Metadata File (.txt/.tsv/.csv)", accept = c(".txt", ".tsv", ".csv")),
      helpText("Metadata must contain columns: SampleID and Group"),
      hr(),

      h4("2. Plot Controls"),
      sliderInput("plot_height", "Global Plot Height (px):", min = 400, max = 3000, value = 800),
      sliderInput("plot_width", "Global Plot Width (%):", min = 50, max = 100, value = 100),
      sliderInput("top_n", "Top N taxa for plots:", min = 5, max = 300, value = 50, step = 5),
      checkboxInput("show_all_sig_species", "Show all significant species/biomarkers", value = TRUE),
      checkboxInput("collapse_other_phyla", "Collapse non-top phyla as Other", value = TRUE),
      checkboxInput("log_heatmap", "Use log10(x + 1) scale for heatmaps", value = FALSE),
      radioButtons(
        "heatmap_scale",
        "Species / ANI heatmap abundance scale:",
        choices = c(
          "Normalized Percent Abundance" = "percent",
          "Normalized Counts" = "norm_counts"
        ),
        selected = "percent"
      ),
      selectInput(
        "ordination_method",
        "Ordination method:",
        choices = c(
          "Bray-Curtis PCoA" = "bray_pcoa",
          "Hellinger PCA" = "hellinger_pca"
        ),
        selected = "bray_pcoa"
      ),
      hr(),

      h4("3. Biomarker Filters"),
      sliderInput("min_prevalence", "Minimum prevalence (% samples):", min = 0, max = 100, value = 10, step = 5),
      numericInput("min_mean_count", "Minimum mean normalized count:", value = 0, min = 0, step = 1),
      numericInput("lda_cutoff", "Minimum LDA-like effect size:", value = 0, min = 0, step = 0.1),
      hr(),

      h4("4. Analysis"),
      actionButton("run", "Run Analysis", class = "btn-success btn-lg", width = "100%"),
      hr(),

      h4("5. Export"),
      downloadButton("download_report", "Download HTML Plot Report", class = "btn-primary", width = "100%"),
      br(), br(),
      downloadButton("download_norm_reads", "Download Normalized Reads (.txt)", class = "btn-info", width = "100%"),
      br(), br(),
      downloadButton("download_percent", "Download Percent Abundance (.txt)", class = "btn-info", width = "100%"),
      br(), br(),
      downloadButton("download_biomarkers", "Download Biomarker Table (.txt)", class = "btn-warning", width = "100%"),
      br(), br(),
      downloadButton("download_ani100", "Download ANI 100 Table (.txt)", class = "btn-secondary", width = "100%"),
      br(), br(),
      downloadButton("download_session_info", "Download Session Info (.txt)", class = "btn-secondary", width = "100%"),
      br(), br(),
      downloadButton("download_citation", "Download CITATION.cff", class = "btn-secondary", width = "100%"),
      br(), br(),
      selectInput(
        "plot_to_export",
        "Static plot export:",
        choices = c(
          "Phylum composition" = "phylum",
          "Species heatmap" = "species_heatmap",
          "Ordination" = "ordination",
          "Bray-Curtis dendrogram" = "dendrogram",
          "Phylum biomarkers" = "biomarker_phylum",
          "Species biomarkers" = "biomarker_species",
          "ANI 100 heatmap" = "ani100"
        )
      ),
      numericInput("export_width", "Export width (in):", value = 12, min = 4, max = 40, step = 1),
      numericInput("export_height", "Export height (in):", value = 8, min = 4, max = 60, step = 1),
      numericInput("export_dpi", "PNG resolution (dpi):", value = 300, min = 72, max = 600, step = 50),
      downloadButton("download_selected_plot_png", "Download Selected Plot (.png)", class = "btn-success", width = "100%"),
      br(), br(),
      downloadButton("download_selected_plot_pdf", "Download Selected Plot (.pdf)", class = "btn-success", width = "100%"),
      hr(),

      radioButtons(
        "tax_level", "Taxonomic Rank View:",
        choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
        selected = "Phylum"
      )
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Grouping", br(), DTOutput("group_table")),

        tabPanel(
          "Clean Taxonomy Table", br(),
          helpText("Normalized percentage abundances aggregated by selected taxonomic rank."),
          DTOutput("dynamic_table")
        ),

        tabPanel(
          "Composition", br(),
          h4("Phylum Composition from Normalized Data"),
          uiOutput("phylum_ui"), hr(),
          h4("Top Species Heatmap"),
          uiOutput("species_ui")
        ),

        tabPanel(
          "Beta Diversity", br(),
          uiOutput("pca_ui"), hr(),
          h4("PERMANOVA"),
          DTOutput("permanova_table"),
          hr(),
          h4("Interactive Bray-Curtis Dendrogram"),
          uiOutput("dendro_ui")
        ),

        tabPanel(
          "Biomarkers", br(),
          helpText("Biomarker discovery uses normalized read counts with Kruskal-Wallis screening, one-vs-rest Wilcoxon validation, and LDA-like effect-size estimation. This is LEfSe-like, not the original LEfSe implementation."),
          h4("Phylum-level biomarkers"),
          uiOutput("lefse_phylum_ui"),
          hr(),
          h4("Species-level biomarkers"),
          uiOutput("lefse_species_ui"),
          hr(),
          h4("Biomarker result table"),
          DTOutput("biomarker_table")
        ),

        tabPanel(
          "ANI 100", br(),
          helpText("Species and strain-level matches with ANI = 100."),
          h4("ANI 100 abundance heatmap"),
          uiOutput("ani100_ui"),
          hr(),
          h4("ANI 100 table"),
          DTOutput("ani100_table")
        ),

        tabPanel(
          "About / Reproducibility", br(),
          h4("About this app"),
          p(strong("Sylph MPA Metagenomics Explorer"), " is an R Shiny application for interactive analysis of Sylph MPA taxonomic profiles."),
          tags$ul(
            tags$li("Normalises sequence abundance across samples using sample-depth scaling."),
            tags$li("Retains MPA-style strain filtering to reduce species/strain double counting."),
            tags$li("Provides composition plots, heatmaps, ordination, PERMANOVA, LEfSe-like biomarkers, and ANI = 100 summaries."),
            tags$li("Exports processed tables, plots, reports, citation metadata, and session information.")
          ),
          h4("Recommended citation"),
          verbatimTextOutput("citation_text"),
          h4("Session information"),
          verbatimTextOutput("session_info")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  v <- reactiveValues(
    raw_long_all = NULL,
    count_df = NULL,
    percent_df = NULL,
    group_df = NULL,
    ps = NULL,
    ps_counts = NULL,
    biomarker_table = NULL,
    ani100_table = NULL,
    permanova_table = NULL,
    n_samples = NULL,
    validation_messages = NULL
  )

  clean_tax <- function(x) {
    x <- stringr::str_remove(x, "^[dkpcofgst]__")
    x <- ifelse(is.na(x) | x == "" | x == "NA", "Unclassified", x)
    x
  }

  group_palette <- function(groups) {
    groups <- sort(unique(as.character(groups)))
    cols <- grDevices::hcl.colors(max(length(groups), 1), palette = "Dark 3")
    setNames(cols[seq_along(groups)], groups)
  }

  taxon_palette <- function(n) {
    n <- max(n, 1)
    if (n <= 22) {
      as.character(pals::alphabet(n))
    } else if (n <= 34) {
      as.character(pals::alphabet2(n))
    } else {
      rep(as.character(pals::stepped(48)), length.out = n)
    }
  }

  parse_sylph_advanced <- function(path, fname) {
    header <- readLines(path, n = 1, warn = FALSE)
    sample_id <- stringr::str_extract(header, "(?<=#SampleID\\t)[^\\t]+")
    if (is.na(sample_id) || sample_id == "") {
      sample_id <- gsub("\\.sylphmpa$|\\.tsv$|\\.txt$|\\.csv$", "", fname, ignore.case = TRUE)
    }

    dat <- readr::read_tsv(path, comment = "#", show_col_types = FALSE, progress = FALSE)
    if (nrow(dat) == 0) stop("The file contains no data rows after header/comment removal: ", fname)
    cn <- colnames(dat)

    clade_col <- cn[stringr::str_detect(cn, regex("^clade_name$", ignore_case = TRUE))][1]
    seq_col   <- cn[stringr::str_detect(cn, regex("^sequence_abundance$", ignore_case = TRUE))][1]
    ani_col   <- cn[stringr::str_detect(cn, regex("^ANI$|ANI", ignore_case = TRUE))][1]
    cov_col   <- cn[stringr::str_detect(cn, regex("^Coverage$|Coverage", ignore_case = TRUE))][1]

    if (is.na(clade_col)) stop("No clade_name column found in: ", fname)
    if (is.na(seq_col)) stop("No sequence_abundance column found in: ", fname)

    out <- dat %>%
      transmute(
        SampleID = sample_id,
        clade = as.character(.data[[clade_col]]),
        sequence_abundance = suppressWarnings(as.numeric(.data[[seq_col]])),
        ANI = if (!is.na(ani_col)) suppressWarnings(as.numeric(.data[[ani_col]])) else NA_real_,
        Coverage = if (!is.na(cov_col)) suppressWarnings(as.numeric(.data[[cov_col]])) else NA_real_
      ) %>%
      filter(!is.na(clade), clade != "", !is.na(sequence_abundance))

    if (nrow(out) == 0) stop("No valid clade/sequence_abundance rows found in: ", fname)
    if (any(out$sequence_abundance < 0, na.rm = TRUE)) stop("Negative sequence_abundance values found in: ", fname)
    out
  }

  split_taxonomy <- function(clade_vec) {
    tax_mat <- stringr::str_split_fixed(clade_vec, "\\|", 7) %>%
      as.data.frame(stringsAsFactors = FALSE)
    colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    tax_mat %>% mutate(across(everything(), clean_tax))
  }

  build_phyloseq_from_wide <- function(wide_df, group_df) {
    sp_df <- wide_df %>%
      filter(stringr::str_detect(clade, "s__")) %>%
      filter(!stringr::str_detect(clade, "t__"))

    validate(need(nrow(sp_df) > 0, "No species-level rows without strain labels were found."))

    otu_mat <- as.matrix(sp_df[, -1, drop = FALSE])
    rownames(otu_mat) <- sp_df$clade
    storage.mode(otu_mat) <- "numeric"

    tax_mat <- split_taxonomy(sp_df$clade)
    rownames(tax_mat) <- sp_df$clade

    meta_df <- group_df %>% distinct(SampleID, .keep_all = TRUE)
    rownames(meta_df) <- meta_df$SampleID

    common_samples <- intersect(colnames(otu_mat), meta_df$SampleID)
    validate(need(length(common_samples) >= 1, "No matching SampleID values between abundance table and metadata."))

    otu_mat <- otu_mat[, common_samples, drop = FALSE]
    meta_df <- meta_df[common_samples, , drop = FALSE]

    phyloseq(
      otu_table(otu_mat, taxa_are_rows = TRUE),
      tax_table(as.matrix(tax_mat)),
      sample_data(meta_df)
    )
  }

  make_species_tax_table <- function(wide_df) {
    tax_df <- split_taxonomy(wide_df$clade)
    bind_cols(tax_df, wide_df %>% dplyr::select(-clade))
  }


  app_citation_text <- function() {
    paste0(
      "Wirth R. SylphScope: an interactive Shiny application ",
      "for exploratory analysis of Sylph MPA taxonomic profiles. Version ", APP_VERSION,
      ". Available from GitHub. Manuscript in preparation."
    )
  }

  citation_cff_text <- function() {
    paste(
      "cff-version: 1.2.0",
      "message: If you use this software, please cite it using the metadata from this file.",
      "title: Sylph MPA Metagenomics Explorer",
      paste0("version: ", APP_VERSION),
      "type: software",
      "authors:",
      "  - family-names: Wirth",
      "    given-names: Roland",
      "abstract: Interactive R Shiny application for exploratory analysis, visualisation and statistical comparison of Sylph MPA taxonomic abundance tables.",
      "keywords:",
      "  - metagenomics",
      "  - microbiome",
      "  - Sylph",
      "  - Shiny",
      "  - microbial ecology",
      sep = "\n"
    )
  }

  validate_uploaded_dataset <- function(raw_long_all) {
    msgs <- character()
    if (anyDuplicated(raw_long_all %>% dplyr::select(SampleID, clade)) > 0) {
      msgs <- c(msgs, "Duplicate SampleID/clade rows were detected and will be summed during matrix construction.")
    }
    if (any(!stringr::str_detect(raw_long_all$clade, "[dkpcofgst]__"))) {
      msgs <- c(msgs, "Some clade names do not contain standard MPA prefixes such as p__, g__, s__, or t__.")
    }
    if (length(unique(raw_long_all$SampleID)) == 1) {
      msgs <- c(msgs, "Single-sample mode: beta diversity, PERMANOVA and group biomarkers require multiple samples/groups.")
    }
    if (length(msgs) == 0) msgs <- "No major input problems detected."
    msgs
  }

  safe_ggsave <- function(plot_obj, file, width, height, dpi, device = NULL) {
    ggplot2::ggsave(
      filename = file, plot = plot_obj, width = width, height = height,
      dpi = dpi, units = "in", limitsize = FALSE, device = device
    )
  }

  filter_prevalent_taxa <- function(ps_obj, min_prev_percent = 0, min_mean_count = 0) {
    if (ntaxa(ps_obj) == 0) return(ps_obj)
    otu <- as(otu_table(ps_obj), "matrix")
    if (!taxa_are_rows(ps_obj)) otu <- t(otu)

    prevalence <- rowSums(otu > 0, na.rm = TRUE) / ncol(otu) * 100
    mean_abund <- rowMeans(otu, na.rm = TRUE)
    keep_taxa <- names(prevalence)[prevalence >= min_prev_percent & mean_abund >= min_mean_count]

    if (length(keep_taxa) == 0) {
      keep_taxa <- names(sort(taxa_sums(ps_obj), decreasing = TRUE))[1:min(10, ntaxa(ps_obj))]
    }
    prune_taxa(keep_taxa, ps_obj)
  }

  make_lefse_like_table <- function(ps_counts, rank = "Species", min_prev = 0, min_mean = 0, lda_cutoff = 0) {
    ps_rank <- tax_glom(ps_counts, rank, NArm = FALSE)
    ps_rank <- filter_prevalent_taxa(ps_rank, min_prev, min_mean)

    dat <- psmelt(ps_rank) %>%
      as_tibble() %>%
      mutate(Group = as.factor(Group), TaxonName = .data[[rank]]) %>%
      filter(!is.na(TaxonName), TaxonName != "", TaxonName != "Unclassified", !is.na(Group))

    if (n_distinct(dat$Group) < 2 || nrow(dat) == 0) return(tibble())

    group_means <- dat %>%
      group_by(OTU, TaxonName, Group) %>%
      summarise(mean_counts = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Group, values_from = mean_counts, names_prefix = "mean_NormCounts_", values_fill = 0)

    kw_table <- dat %>%
      group_by(OTU, TaxonName) %>%
      filter(n_distinct(Group) >= 2) %>%
      summarise(
        kw_pvalue = tryCatch(kruskal.test(Abundance ~ Group)$p.value, error = function(e) NA_real_),
        TopGroup = names(which.max(tapply(Abundance, Group, mean, na.rm = TRUE))),
        MaxMeanCounts = max(tapply(Abundance, Group, mean, na.rm = TRUE)),
        PrevalencePercent = mean(Abundance > 0, na.rm = TRUE) * 100,
        Rank = rank,
        .groups = "drop"
      ) %>%
      filter(!is.na(kw_pvalue)) %>%
      mutate(kw_padj = p.adjust(kw_pvalue, method = "BH"))

    wilcox_lda <- dat %>%
      group_by(OTU, TaxonName) %>%
      group_modify(function(.x, .y) {
        top_group <- names(which.max(tapply(.x$Abundance, .x$Group, mean, na.rm = TRUE)))
        tmp <- .x %>%
          mutate(
            is_top = ifelse(Group == top_group, "TopGroup", "Other"),
            value_log = log10(Abundance + 1)
          )

        wilcox_p <- tryCatch(
          wilcox.test(Abundance ~ is_top, data = tmp, exact = FALSE)$p.value,
          error = function(e) NA_real_
        )

        lda_score <- tryCatch({
          if (n_distinct(tmp$Group) < 2 || sd(tmp$value_log, na.rm = TRUE) == 0) {
            NA_real_
          } else {
            fit <- MASS::lda(Group ~ value_log, data = tmp)
            pred <- predict(fit)$x[, 1]
            abs(mean(pred[tmp$Group == top_group], na.rm = TRUE) - mean(pred[tmp$Group != top_group], na.rm = TRUE))
          }
        }, error = function(e) NA_real_)

        tibble(wilcox_pvalue = wilcox_p, lda_score = lda_score)
      }) %>%
      ungroup()

    kw_table %>%
      left_join(wilcox_lda, by = c("OTU", "TaxonName")) %>%
      mutate(
        wilcox_padj = p.adjust(wilcox_pvalue, method = "BH"),
        Significant = kw_padj < 0.05 & wilcox_padj < 0.05 & !is.na(lda_score) & lda_score >= lda_cutoff
      ) %>%
      left_join(group_means, by = c("OTU", "TaxonName")) %>%
      arrange(desc(Significant), desc(lda_score), kw_padj, wilcox_padj)
  }

  make_lefse_plot_gg <- function(biomarker_table, rank = "Species", top_n = 50, show_all = FALSE) {
    dat <- biomarker_table %>% filter(Rank == rank, Significant, !is.na(lda_score))
    if (!show_all) dat <- dat %>% slice_max(lda_score, n = top_n)
    dat <- dat %>% arrange(lda_score) %>% mutate(TaxonName = factor(TaxonName, levels = unique(TaxonName)))

    if (nrow(dat) == 0) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("No significant", rank, "biomarkers")) +
          theme_void()
      )
    }

    ggplot(dat, aes(
      x = TaxonName,
      y = lda_score,
      fill = TopGroup,
      text = paste0(
        "Taxon: ", TaxonName,
        "<br>Top group: ", TopGroup,
        "<br>LDA-like score: ", signif(lda_score, 3),
        "<br>KW FDR: ", signif(kw_padj, 3),
        "<br>Wilcoxon FDR: ", signif(wilcox_padj, 3)
      )
    )) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      scale_fill_manual(values = group_palette(dat$TopGroup)) +
      labs(title = paste(rank, "biomarkers"), x = NULL, y = "LDA-like effect size", fill = "Enriched group")
  }

  make_static_plot <- function(plot_id) {
    req(v$ps)

    if (plot_id == "phylum") {
      ps_phylum <- tax_glom(v$ps, "Phylum", NArm = FALSE)
      top_taxa <- names(sort(taxa_sums(ps_phylum), decreasing = TRUE))[1:min(input$top_n, ntaxa(ps_phylum))]
      df <- psmelt(ps_phylum) %>% as_tibble()
      if (input$collapse_other_phyla) {
        top_phylum_names <- df %>% filter(OTU %in% top_taxa) %>% distinct(Phylum) %>% pull(Phylum)
        df <- df %>% mutate(Phylum = ifelse(Phylum %in% top_phylum_names, Phylum, "Other")) %>%
          group_by(Sample, Group, Phylum) %>% summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
      } else {
        df <- df %>% filter(OTU %in% top_taxa)
      }
      return(
        ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
          geom_col() + facet_grid(~Group, scales = "free_x", space = "free") +
          scale_fill_manual(values = taxon_palette(n_distinct(df$Phylum))) + theme_bw() +
          labs(title = "Phylum relative abundance from normalised profiles", x = NULL, y = "Relative abundance (%)") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      )
    }

    if (plot_id == "species_heatmap") {
      ps_source <- if (input$heatmap_scale == "percent") v$ps else v$ps_counts
      req(ps_source)
      top_taxa <- names(sort(taxa_sums(ps_source), decreasing = TRUE))[1:min(input$top_n, ntaxa(ps_source))]
      df <- psmelt(prune_taxa(top_taxa, ps_source)) %>% as_tibble() %>%
        mutate(PlotValue = if (input$log_heatmap) log10(Abundance + 1) else Abundance, Species = fct_reorder(Species, Abundance, .fun = sum, .desc = FALSE))
      return(
        ggplot(df, aes(x = Sample, y = Species, fill = PlotValue)) + geom_tile() +
          facet_grid(~Group, scales = "free_x", space = "free") +
          scale_fill_viridis_c(option = "inferno", direction = -1) + theme_bw() +
          labs(title = "Species abundance heatmap", x = NULL, y = NULL, fill = "Abundance") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
      )
    }

    if (plot_id == "ordination") {
      if (input$ordination_method == "bray_pcoa") {
        ord <- ordinate(v$ps, method = "PCoA", distance = "bray")
        df_pca <- data.frame(ord$vectors[, 1:2, drop = FALSE], SampleID = rownames(ord$vectors), Group = sample_data(v$ps)$Group)
        colnames(df_pca)[1:2] <- c("Axis.1", "Axis.2")
        title <- "Bray-Curtis PCoA"; xlab <- "Axis 1"; ylab <- "Axis 2"
      } else {
        otu <- as(otu_table(v$ps), "matrix"); if (taxa_are_rows(v$ps)) otu <- t(otu)
        pca <- prcomp(vegan::decostand(otu, method = "hellinger"), center = TRUE, scale. = FALSE)
        var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
        df_pca <- data.frame(Axis.1 = pca$x[, 1], Axis.2 = pca$x[, 2], SampleID = rownames(pca$x), Group = sample_data(v$ps)$Group)
        title <- "Hellinger PCA"; xlab <- paste0("PC1 (", var_exp[1], "%)"); ylab <- paste0("PC2 (", var_exp[2], "%)")
      }
      return(ggplot(df_pca, aes(Axis.1, Axis.2, color = Group, label = SampleID)) + geom_point(size = 4) + geom_text(vjust = -1.2, size = 3, show.legend = FALSE) + theme_bw() + scale_color_manual(values = group_palette(df_pca$Group)) + labs(title = title, x = xlab, y = ylab))
    }

    if (plot_id == "dendrogram") {
      dist <- phyloseq::distance(v$ps, method = "bray")
      hc <- hclust(dist, method = "average")
      dend_data <- dendextend::as.ggdend(as.dendrogram(hc))
      meta <- data.frame(SampleID = rownames(sample_data(v$ps)), Group = sample_data(v$ps)$Group)
      nodes <- dend_data$labels %>% left_join(meta, by = c("label" = "SampleID"))
      return(ggplot() + geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) + geom_point(data = nodes, aes(x = x, y = y, color = Group), size = 3) + scale_color_manual(values = group_palette(nodes$Group)) + theme_void() + labs(title = "Bray-Curtis dendrogram, average linkage"))
    }

    if (plot_id == "biomarker_phylum") { req(v$biomarker_table); return(make_lefse_plot_gg(v$biomarker_table, "Phylum", input$top_n, FALSE)) }
    if (plot_id == "biomarker_species") { req(v$biomarker_table); return(make_lefse_plot_gg(v$biomarker_table, "Species", input$top_n, input$show_all_sig_species)) }

    if (plot_id == "ani100") {
      req(v$ani100_table)
      if (nrow(v$ani100_table) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No ANI = 100 hits found") + theme_void())
      top_labels <- v$ani100_table %>% group_by(Label) %>% summarise(total_abundance = sum(norm_percent, na.rm = TRUE), .groups = "drop") %>% arrange(desc(total_abundance)) %>% slice_head(n = input$top_n) %>% pull(Label)
      df <- v$ani100_table %>% filter(Label %in% top_labels) %>% mutate(PlotAbundance = if (input$heatmap_scale == "percent") norm_percent else norm_reads, PlotValue = if (input$log_heatmap) log10(PlotAbundance + 1) else PlotAbundance, Label = fct_reorder(Label, PlotAbundance, .fun = sum, .desc = FALSE))
      return(ggplot(df, aes(SampleID, Label, fill = PlotValue)) + geom_tile() + scale_fill_viridis_c(option = "inferno", direction = -1) + theme_bw() + labs(title = "ANI = 100 species / strain matches", x = NULL, y = NULL, fill = "Abundance") + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    }

    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Selected plot is not available") + theme_void()
  }

  observeEvent(input$files, {
    req(input$files)

    parsed <- tryCatch({
      purrr::map_dfr(1:nrow(input$files), function(i) {
        parse_sylph_advanced(input$files$datapath[i], input$files$name[i])
      })
    }, error = function(e) {
      showNotification(paste("File parsing failed:", e$message), type = "error", duration = 10)
      NULL
    })

    req(parsed)
    validate(need(nrow(parsed) > 0, "No usable rows found in uploaded Sylph MPA files."))

    raw_long_all <- parsed %>%
      mutate(SampleID = as.character(SampleID))

    # MPA-style strain filtering retained intentionally:
    # if a species has strain-level rows anywhere in the dataset, the species summary row is excluded from the depth total.
    # This prevents double counting of species and strain rows.
    depth_base <- raw_long_all %>%
      filter(str_detect(clade, "s__")) %>%
      mutate(
        has_strain = str_detect(clade, "t__"),
        species_part = str_extract(clade, "s__[^|]+")
      )

    species_with_strains <- depth_base %>%
      filter(has_strain) %>%
      pull(species_part) %>%
      unique()

    clean_reference_rows <- depth_base %>%
      filter(has_strain | !(species_part %in% species_with_strains))

    sample_sums <- clean_reference_rows %>%
      group_by(SampleID) %>%
      summarise(total_sum = sum(sequence_abundance, na.rm = TRUE), .groups = "drop")

    avg_sum <- mean(sample_sums$total_sum[sample_sums$total_sum > 0], na.rm = TRUE)
    validate(need(is.finite(avg_sum) && avg_sum > 0, "Normalized depth could not be calculated. Check sequence_abundance values."))

    normalized_matrix <- raw_long_all %>%
      left_join(sample_sums, by = "SampleID") %>%
      mutate(
        total_sum = replace_na(total_sum, 0),
        norm_reads = ifelse(total_sum > 0, sequence_abundance * (avg_sum / total_sum), 0),
        norm_percent = ifelse(total_sum > 0, (norm_reads / avg_sum) * 100, 0)
      ) %>%
      ungroup()

    v$raw_long_all <- normalized_matrix
    v$n_samples <- n_distinct(normalized_matrix$SampleID)
    v$validation_messages <- validate_uploaded_dataset(normalized_matrix)
    showNotification(paste(v$validation_messages, collapse = " | "), type = "message", duration = 8)

    species_only_matrix <- normalized_matrix %>%
      filter(str_detect(clade, "s__"), !str_detect(clade, "t__"))

    v$count_df <- species_only_matrix %>%
      dplyr::select(clade, SampleID, norm_reads) %>%
      group_by(clade, SampleID) %>%
      summarise(norm_reads = sum(norm_reads, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = SampleID, values_from = norm_reads, values_fill = 0)

    v$percent_df <- species_only_matrix %>%
      dplyr::select(clade, SampleID, norm_percent) %>%
      group_by(clade, SampleID) %>%
      summarise(norm_percent = sum(norm_percent, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = SampleID, values_from = norm_percent, values_fill = 0)

    v$group_df <- data.frame(
      SampleID = unique(normalized_matrix$SampleID),
      Group = "Group_1",
      stringsAsFactors = FALSE
    )

    v$ps <- NULL
    v$ps_counts <- NULL
    v$biomarker_table <- NULL
    v$ani100_table <- NULL
    v$permanova_table <- NULL
  })

  observeEvent(input$metadata_file, {
    req(input$metadata_file, v$group_df)
    ext <- tolower(tools::file_ext(input$metadata_file$name))

    meta <- tryCatch({
      if (ext == "csv") readr::read_csv(input$metadata_file$datapath, show_col_types = FALSE)
      else readr::read_tsv(input$metadata_file$datapath, show_col_types = FALSE)
    }, error = function(e) {
      showNotification(paste("Metadata reading failed:", e$message), type = "error", duration = 10)
      NULL
    })

    req(meta)
    validate(need(all(c("SampleID", "Group") %in% colnames(meta)), "Metadata must contain SampleID and Group."))

    v$group_df <- v$group_df %>%
      dplyr::select(SampleID) %>%
      left_join(meta %>% dplyr::select(SampleID, Group) %>% distinct(SampleID, .keep_all = TRUE), by = "SampleID") %>%
      mutate(Group = ifelse(is.na(Group), "Group_1", as.character(Group)))
  })

  output$group_table <- renderDT({
    req(v$group_df)
    datatable(
      v$group_df,
      rownames = FALSE,
      editable = list(target = "cell", disable = list(columns = c(0))),
      options = list(dom = "t", pageLength = 100)
    )
  })

  observeEvent(input$group_table_cell_edit, {
    info <- input$group_table_cell_edit
    if (info$col + 1 == which(colnames(v$group_df) == "Group")) {
      v$group_df[info$row, "Group"] <- as.character(info$value)
    }
  })

  observeEvent(input$run, {
    req(v$percent_df, v$count_df, v$group_df)

    withProgress(message = "Processing calculations...", {
      incProgress(0.2, detail = "Building phyloseq objects")
      v$ps <- build_phyloseq_from_wide(v$percent_df, v$group_df)
      v$ps_counts <- build_phyloseq_from_wide(v$count_df, v$group_df)

      incProgress(0.25, detail = "Running biomarker analysis")
      if (n_distinct(v$group_df$Group) >= 2 && nrow(v$group_df) >= 3) {
        v$biomarker_table <- bind_rows(
          make_lefse_like_table(v$ps_counts, "Phylum", input$min_prevalence, input$min_mean_count, input$lda_cutoff),
          make_lefse_like_table(v$ps_counts, "Species", input$min_prevalence, input$min_mean_count, input$lda_cutoff)
        )
      } else {
        v$biomarker_table <- tibble()
      }

      incProgress(0.25, detail = "Extracting ANI 100 entries")
      v$ani100_table <- v$raw_long_all %>%
        filter(!is.na(ANI), ANI == 100) %>%
        mutate(
          Species = clean_tax(str_extract(clade, "s__[^|]+")),
          Strain = clean_tax(str_extract(clade, "t__[^|]+")),
          Strain = ifelse(Strain == "Unclassified", NA_character_, Strain),
          Label = case_when(!is.na(Strain) ~ paste0(Species, " | ", Strain), TRUE ~ Species)
        ) %>%
        dplyr::select(SampleID, Species, Strain, Label, ANI, Coverage, norm_reads, norm_percent, clade)

      incProgress(0.2, detail = "Running PERMANOVA")
      if (n_distinct(sample_data(v$ps)$Group) >= 2 && nsamples(v$ps) >= 3) {
        dist <- phyloseq::distance(v$ps, method = "bray")
        meta <- data.frame(sample_data(v$ps))
        v$permanova_table <- tryCatch({
          as.data.frame(vegan::adonis2(dist ~ Group, data = meta, permutations = 999)) %>%
            rownames_to_column("Term")
        }, error = function(e) tibble(Term = "PERMANOVA failed", Message = e$message))
      } else {
        v$permanova_table <- tibble(Message = "PERMANOVA requires at least two groups and at least three samples.")
      }
    })
  })

  output$phylum_ui       <- renderUI({ plotlyOutput("phylum_bar", height = input$plot_height, width = paste0(input$plot_width, "%")) })
  output$species_ui      <- renderUI({ plotlyOutput("species_heatmap", height = max(input$plot_height, input$top_n * 22), width = paste0(input$plot_width, "%")) })
  output$pca_ui          <- renderUI({ plotlyOutput("pca_plot", height = input$plot_height, width = paste0(input$plot_width, "%")) })
  output$dendro_ui       <- renderUI({ plotlyOutput("dendrogram", height = input$plot_height, width = paste0(input$plot_width, "%")) })
  output$lefse_phylum_ui <- renderUI({ plotlyOutput("lefse_phylum_plot", height = input$plot_height, width = paste0(input$plot_width, "%")) })

  output$lefse_species_ui <- renderUI({
    req(v$biomarker_table)
    n_sig <- v$biomarker_table %>% filter(Rank == "Species", Significant) %>% nrow()
    if (!input$show_all_sig_species) n_sig <- min(n_sig, input$top_n)
    plotlyOutput("lefse_species_plot", height = max(input$plot_height, n_sig * 28), width = paste0(input$plot_width, "%"))
  })

  output$ani100_ui <- renderUI({
    req(v$ani100_table)
    plotlyOutput("ani100_plot", height = max(input$plot_height, min(n_distinct(v$ani100_table$Label), input$top_n) * 24), width = paste0(input$plot_width, "%"))
  })

  output$phylum_bar <- renderPlotly({
    req(v$ps)
    ps_phylum <- tax_glom(v$ps, "Phylum", NArm = FALSE)
    top_taxa <- names(sort(taxa_sums(ps_phylum), decreasing = TRUE))[1:min(input$top_n, ntaxa(ps_phylum))]

    df <- psmelt(ps_phylum) %>% as_tibble()

    if (input$collapse_other_phyla) {
      top_phylum_names <- df %>%
        filter(OTU %in% top_taxa) %>%
        distinct(Phylum) %>%
        pull(Phylum)

      df <- df %>%
        mutate(Phylum = ifelse(Phylum %in% top_phylum_names, Phylum, "Other")) %>%
        group_by(Sample, Group, Phylum) %>%
        summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
    } else {
      df <- df %>% filter(OTU %in% top_taxa)
    }

    df <- df %>% mutate(Phylum = factor(Phylum, levels = unique(Phylum[order(-Abundance)])))

    p <- ggplot(df, aes(
      x = Sample,
      y = Abundance,
      fill = Phylum,
      text = paste0("Sample: ", Sample, "<br>Phylum: ", Phylum, "<br>Abundance: ", signif(Abundance, 4), "%")
    )) +
      geom_col(position = "stack") +
      facet_grid(~Group, scales = "free_x", space = "free") +
      scale_fill_manual(values = taxon_palette(n_distinct(df$Phylum))) +
      theme_bw() +
      labs(title = "Phylum Relative Abundance from Normalized Profiles", x = NULL, y = "Relative Abundance (%)") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggplotly(p, tooltip = "text")
  })

  output$species_heatmap <- renderPlotly({
    req(v$ps, v$ps_counts)

    if (input$heatmap_scale == "percent") {
      top_taxa <- names(sort(taxa_sums(v$ps), decreasing = TRUE))[1:min(input$top_n, ntaxa(v$ps))]
      df <- psmelt(prune_taxa(top_taxa, v$ps)) %>% mutate(PlotAbundance = Abundance)
      fill_lab <- "Normalized Percent (%)"
    } else {
      top_taxa <- names(sort(taxa_sums(v$ps_counts), decreasing = TRUE))[1:min(input$top_n, ntaxa(v$ps_counts))]
      df <- psmelt(prune_taxa(top_taxa, v$ps_counts)) %>% mutate(PlotAbundance = Abundance)
      fill_lab <- "Normalized Counts"
    }

    df <- df %>%
      as_tibble() %>%
      mutate(
        PlotValue = if (input$log_heatmap) log10(PlotAbundance + 1) else PlotAbundance,
        Species = fct_reorder(Species, PlotAbundance, .fun = sum, .desc = FALSE)
      )

    if (input$log_heatmap) fill_lab <- paste0("log10(", fill_lab, " + 1)")

    p <- ggplot(df, aes(
      x = Sample,
      y = Species,
      fill = PlotValue,
      text = paste0("Sample: ", Sample, "<br>Species: ", Species, "<br>Value: ", signif(PlotAbundance, 4))
    )) +
      geom_tile() +
      facet_grid(~Group, scales = "free_x", space = "free") +
      scale_fill_viridis_c(option = "inferno", direction = -1) +
      theme_bw() +
      labs(title = "Species Abundance Heatmap", x = NULL, y = NULL, fill = fill_lab) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggplotly(p, tooltip = "text")
  })

  output$pca_plot <- renderPlotly({
    req(v$ps)

    if (input$ordination_method == "bray_pcoa") {
      ord <- ordinate(v$ps, method = "PCoA", distance = "bray")
      df_pca <- data.frame(ord$vectors[, 1:2, drop = FALSE], SampleID = rownames(ord$vectors), Group = sample_data(v$ps)$Group)
      colnames(df_pca)[1:2] <- c("Axis.1", "Axis.2")
      title <- "Bray-Curtis PCoA"
      xlab <- "Axis 1"
      ylab <- "Axis 2"
    } else {
      otu <- as(otu_table(v$ps), "matrix")
      if (taxa_are_rows(v$ps)) otu <- t(otu)
      otu_hel <- vegan::decostand(otu, method = "hellinger")
      pca <- prcomp(otu_hel, center = TRUE, scale. = FALSE)
      var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
      df_pca <- data.frame(Axis.1 = pca$x[, 1], Axis.2 = pca$x[, 2], SampleID = rownames(pca$x), Group = sample_data(v$ps)$Group)
      title <- "Hellinger PCA"
      xlab <- paste0("PC1 (", var_exp[1], "%)")
      ylab <- paste0("PC2 (", var_exp[2], "%)")
    }

    p <- ggplot(df_pca, aes(x = Axis.1, y = Axis.2, color = Group, text = SampleID)) +
      geom_point(size = 4) +
      geom_text(aes(label = SampleID), vjust = -1.2, size = 3, show.legend = FALSE) +
      theme_bw() +
      scale_color_manual(values = group_palette(df_pca$Group)) +
      labs(title = title, x = xlab, y = ylab)

    ggplotly(p, tooltip = "text")
  })

  output$permanova_table <- renderDT({
    req(v$permanova_table)
    datatable(v$permanova_table, rownames = FALSE, options = list(scrollX = TRUE, dom = "t", pageLength = 10))
  })

  output$dendrogram <- renderPlotly({
    req(v$ps)
    dist <- phyloseq::distance(v$ps, method = "bray")
    hc <- hclust(dist, method = "average")
    dend_data <- dendextend::as.ggdend(as.dendrogram(hc))
    meta <- data.frame(SampleID = rownames(sample_data(v$ps)), Group = sample_data(v$ps)$Group)
    nodes <- dend_data$labels %>% left_join(meta, by = c("label" = "SampleID"))

    p <- ggplot() +
      geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_point(data = nodes, aes(x = x, y = y, color = Group, text = label), size = 3) +
      scale_color_manual(values = group_palette(nodes$Group)) +
      theme_void() +
      labs(title = "Bray-Curtis dendrogram, average linkage")

    ggplotly(p, tooltip = "text")
  })

  output$lefse_phylum_plot <- renderPlotly({
    req(v$biomarker_table)
    ggplotly(make_lefse_plot_gg(v$biomarker_table, rank = "Phylum", top_n = input$top_n, show_all = FALSE), tooltip = "text")
  })

  output$lefse_species_plot <- renderPlotly({
    req(v$biomarker_table)
    ggplotly(make_lefse_plot_gg(v$biomarker_table, rank = "Species", top_n = input$top_n, show_all = input$show_all_sig_species), tooltip = "text")
  })

  output$biomarker_table <- renderDT({
    req(v$biomarker_table)
    datatable(
      v$biomarker_table %>% filter(Significant),
      extensions = "Buttons",
      rownames = FALSE,
      options = list(scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel"), pageLength = 25)
    )
  })

  output$ani100_table <- renderDT({
    req(v$ani100_table)
    datatable(
      v$ani100_table,
      extensions = "Buttons",
      rownames = FALSE,
      options = list(scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel"), pageLength = 25)
    )
  })

  output$ani100_plot <- renderPlotly({
    req(v$ani100_table)

    if (nrow(v$ani100_table) == 0) {
      return(ggplotly(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No ANI = 100 hits found") + theme_void()))
    }

    top_labels <- v$ani100_table %>%
      group_by(Label) %>%
      summarise(total_abundance = sum(norm_percent, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total_abundance)) %>%
      slice_head(n = input$top_n) %>%
      pull(Label)

    df <- v$ani100_table %>%
      filter(Label %in% top_labels) %>%
      mutate(
        PlotAbundance = if (input$heatmap_scale == "percent") norm_percent else norm_reads,
        PlotValue = if (input$log_heatmap) log10(PlotAbundance + 1) else PlotAbundance,
        Label = fct_reorder(Label, PlotAbundance, .fun = sum, .desc = FALSE)
      )

    fill_lab <- ifelse(input$heatmap_scale == "percent", "Normalized Percent (%)", "Normalized Counts")
    if (input$log_heatmap) fill_lab <- paste0("log10(", fill_lab, " + 1)")

    p <- ggplot(df, aes(
      x = SampleID,
      y = Label,
      fill = PlotValue,
      text = paste0("Sample: ", SampleID, "<br>Label: ", Label, "<br>Value: ", signif(PlotAbundance, 4), "<br>Coverage: ", signif(Coverage, 4))
    )) +
      geom_tile() +
      scale_fill_viridis_c(option = "inferno", direction = -1) +
      theme_bw() +
      labs(title = "ANI = 100 species / strain matches", x = NULL, y = NULL, fill = fill_lab) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggplotly(p, tooltip = "text")
  })

  output$dynamic_table <- renderDT({
    req(v$ps)
    tax_df <- as.data.frame(tax_table(v$ps))
    otu_df <- as.data.frame(otu_table(v$ps))

    if (!taxa_are_rows(v$ps)) otu_df <- as.data.frame(t(otu_df))

    res <- bind_cols(tax_df, otu_df) %>%
      as_tibble() %>%
      group_by(.data[[input$tax_level]]) %>%
      summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
      mutate(TotalAbundance = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
      arrange(desc(TotalAbundance))

    datatable(res, extensions = "Buttons", rownames = FALSE, options = list(scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel"), pageLength = 25))
  })

  output$download_norm_reads <- downloadHandler(
    filename = function() paste0("Species_level_Normalized_reads_", Sys.Date(), ".txt"),
    content = function(file) { req(v$count_df); make_species_tax_table(v$count_df) %>% readr::write_tsv(file) }
  )

  output$download_percent <- downloadHandler(
    filename = function() paste0("Species_level_percent_abundance_", Sys.Date(), ".txt"),
    content = function(file) { req(v$percent_df); make_species_tax_table(v$percent_df) %>% readr::write_tsv(file) }
  )

  output$download_biomarkers <- downloadHandler(
    filename = function() paste0("LEfSe_like_biomarker_table_", Sys.Date(), ".txt"),
    content = function(file) { req(v$biomarker_table); v$biomarker_table %>% filter(Significant) %>% readr::write_tsv(file) }
  )

  output$download_ani100 <- downloadHandler(
    filename = function() paste0("ANI100_species_strain_table_", Sys.Date(), ".txt"),
    content = function(file) { req(v$ani100_table); v$ani100_table %>% readr::write_tsv(file) }
  )


  output$citation_text <- renderText({ app_citation_text() })

  output$session_info <- renderPrint({
    cat(APP_NAME, "\n")
    cat("Version:", APP_VERSION, "\n")
    cat("Date:", as.character(Sys.Date()), "\n\n")
    if (!is.null(v$validation_messages)) {
      cat("Input validation messages:\n")
      cat(paste0("- ", v$validation_messages, collapse = "\n"), "\n\n")
    }
    sessionInfo()
  })

  output$download_session_info <- downloadHandler(
    filename = function() paste0("Sylph_MPA_Explorer_session_info_", Sys.Date(), ".txt"),
    content = function(file) {
      sink(file)
      cat(APP_NAME, "\n")
      cat("Version:", APP_VERSION, "\n")
      cat("Date:", as.character(Sys.Date()), "\n\n")
      if (!is.null(v$validation_messages)) {
        cat("Input validation messages:\n")
        cat(paste0("- ", v$validation_messages, collapse = "\n"), "\n\n")
      }
      print(sessionInfo())
      sink()
    }
  )

  output$download_citation <- downloadHandler(
    filename = function() "CITATION.cff",
    content = function(file) writeLines(citation_cff_text(), file)
  )

  output$download_selected_plot_png <- downloadHandler(
    filename = function() paste0("Sylph_MPA_Explorer_", input$plot_to_export, "_", Sys.Date(), ".png"),
    content = function(file) {
      p <- make_static_plot(input$plot_to_export)
      safe_ggsave(p, file, input$export_width, input$export_height, input$export_dpi, device = "png")
    }
  )

  output$download_selected_plot_pdf <- downloadHandler(
    filename = function() paste0("Sylph_MPA_Explorer_", input$plot_to_export, "_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- make_static_plot(input$plot_to_export)
      safe_ggsave(p, file, input$export_width, input$export_height, input$export_dpi, device = grDevices::cairo_pdf)
    }
  )

  output$download_report <- downloadHandler(
    filename = function() paste0("Sylph_MPA_Explorer_Report_", Sys.Date(), ".html"),
    content = function(file) {
      req(v$ps)
      tmpdir <- tempdir()
      rmd <- file.path(tmpdir, "sylph_report.Rmd")
      saveRDS(v$ps, file.path(tmpdir, "ps.rds"))
      saveRDS(v$biomarker_table, file.path(tmpdir, "biomarkers.rds"))
      saveRDS(v$ani100_table, file.path(tmpdir, "ani100.rds"))
      saveRDS(v$permanova_table, file.path(tmpdir, "permanova.rds"))
      saveRDS(v$validation_messages, file.path(tmpdir, "validation_messages.rds"))

      writeLines(c(
        "---",
        "title: 'Sylph MPA Explorer Report'",
        "output: html_document",
        "---",
        "",
        "```{r setup, message=FALSE, warning=FALSE}",
        "library(tidyverse); library(phyloseq); library(vegan); library(DT)",
        "ps <- readRDS('ps.rds')",
        "biomarkers <- readRDS('biomarkers.rds')",
        "ani100 <- readRDS('ani100.rds')",
        "permanova <- readRDS('permanova.rds')",
        "validation_messages <- readRDS('validation_messages.rds')",
        "```",
        "",
        "## Dataset summary",
        "```{r}",
        "data.frame(Samples = nsamples(ps), Species = ntaxa(ps))",
        "```",
        "",
        "## Input validation messages",
        "```{r}",
        "validation_messages",
        "```",
        "",
        "## PERMANOVA",
        "```{r}",
        "permanova",
        "```",
        "",
        "## Significant biomarkers",
        "```{r}",
        "biomarkers %>% filter(Significant) %>% arrange(Rank, desc(lda_score)) %>% head(100)",
        "```",
        "",
        "## ANI = 100 entries",
        "```{r}",
        "ani100 %>% arrange(desc(norm_percent)) %>% head(100)",
        "```"
      ), rmd)

      oldwd <- setwd(tmpdir)
      on.exit(setwd(oldwd), add = TRUE)
      rmarkdown::render(rmd, output_file = file, quiet = TRUE, envir = new.env(parent = globalenv()))
    }
  )
}

shinyApp(ui, server)
