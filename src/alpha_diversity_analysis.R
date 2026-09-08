################################################################################
# ALPHA DIVERSITY ANALYSIS: FILENAME-BASED LABELING
# Extracts 'year_label' from filename (e.g., marker-LABEL-metric.tsv)
# MODIFIED: Boxplots now feature a Zone-colored background fill with a thin 
# black outline.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(rstatix)
  library(multcompView)
  library(patchwork)
  library(ggsci)      # NPG color palette matching the beta diversity script
})

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

DATA_DIR      <- "E:/projects/Musquash_data/Analysis/results_diversity_reanalysis/no_prvl_filtering/combined_matrices"
METADATA_PATH <- "E:/projects/Musquash_data/Analysis/metadata/comb_station_grp_all_years_metada.tsv"
OUTPUT_DIR    <- "E:/projects/Musquash_data/Analysis/results_diversity_reanalysis/no_prvl_filtering/combined_plots_alpha_diversity"

# Order matches the desired combined-panel order:
analysis_batches <- list(
  list(marker = "mifish", metric_name = "shannon",  label = "Shannon Index"),
  list(marker = "mifish", metric_name = "faith-pd", label = "Faith's PD"),
  list(marker = "COI",    metric_name = "shannon",  label = "Shannon Index"),
  list(marker = "COI",    metric_name = "faith-pd", label = "Faith's PD")
)

# Salinity transparency settings
ALPHA_VAR            <- "not there"
ALPHA_RANGE          <- c(0.10, 0.97)
ALPHA_HIGH_IS_OPAQUE <- TRUE

# ==============================================================================
# 3. HELPER: READ FILE & EXTRACT LABEL FROM NAME
# ==============================================================================

read_labeled_file <- function(filepath) {
  fname      <- basename(filepath)
  parts      <- str_split(fname, "-")[[1]]
  label_val  <- paste0(parts[2], " ",  parts[3])
  
  df <- read_tsv(filepath,
                 col_names = c("sample_id", "diversity_value"),
                 skip = 1, col_types = cols()) %>%
    mutate(year_label = label_val)
  
  return(df)
}

# ==============================================================================
# 4. HELPER: COMPACT LETTER DISPLAY
# ==============================================================================

get_cld_letters <- function(dunn_df) {
  pvals               <- setNames(dunn_df$p.adj,
                                  paste(dunn_df$group1, dunn_df$group2, sep = "-"))
  pvals[is.na(pvals)] <- 1
  tryCatch(
    multcompLetters(pvals)$Letters,
    error = function(e) {
      warning("CLD failed: ", conditionMessage(e))
      NULL
    }
  )
}

# ==============================================================================
# 5. CORE PROCESSING FUNCTION
# ==============================================================================

process_and_plot <- function(data_dir, meta_path, marker, metric_pattern,
                             plot_label, out_dir, save_individual = TRUE) {
  
  message(paste0("\n>>> Processing: ", marker, " - ", plot_label))
  
  # A. FIND & READ FILES
  all_files    <- list.files(data_dir, full.names = TRUE, pattern = "\\.tsv$")
  target_files <- all_files[str_detect(all_files, marker) &
                              str_detect(all_files, metric_pattern)]
  
  if (length(target_files) == 0) return(NULL)
  
  alpha_div <- map_dfr(target_files, read_labeled_file)
  
  # B. LOAD METADATA
  metadata <- read_tsv(meta_path, col_types = cols())
  if (any(str_detect(metadata[1,], "q2:types"))) { metadata <- metadata[-1, ] }
  
  colnames(metadata) <- colnames(metadata) %>% trimws() %>% str_replace_all("-", "_")
  
  if (!"sample_id" %in% colnames(metadata)) {
    id_col_idx <- grep("id", colnames(metadata), ignore.case = TRUE)[1]
    colnames(metadata)[id_col_idx] <- "sample_id"
  }
  
  # C. MERGE & PARSE NUMERIC ALPHA
  custom_order <- c("22 NS_MiU", "23 NS_MiU", "23 MS_MiM", "24 MS_MiM", "25 MS_MiM",
                    "22 NS_LerXT", "23 NS_LerXT", "24 MS_Ler", "25 MS_Ler")
  
  merged_data <- alpha_div %>%
    left_join(metadata, by = "sample_id") %>%
    filter(!is.na(Zone)) %>%
    mutate(
      Zone = as.factor(Zone),
      year_label = factor(year_label, levels = custom_order)
    )
  
  if (nrow(merged_data) == 0) return(NULL)
  
  use_alpha <- ALPHA_VAR %in% colnames(merged_data) && 
    any(!is.na(as.numeric(merged_data[[ALPHA_VAR]])))
  
  if (use_alpha) {
    merged_data[[ALPHA_VAR]] <- as.numeric(merged_data[[ALPHA_VAR]])
  }
  
  # D. STATS & GROUPING
  valid_groups <- merged_data %>%
    group_by(year_label) %>%
    summarise(n_zones = n_distinct(Zone), .groups = "drop") %>%
    filter(n_zones >= 2) %>%
    pull(year_label) %>%
    as.character()
  
  stats_data <- merged_data %>% filter(year_label %in% valid_groups)
  
  kw_results <- stats_data %>%
    group_by(year_label) %>%
    kruskal_test(diversity_value ~ Zone) %>%
    add_significance("p") %>%
    mutate(year_label = factor(year_label, levels = custom_order))
  
  dunn_results <- stats_data %>%
    group_by(year_label) %>%
    dunn_test(diversity_value ~ Zone, p.adjust.method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(year_label = factor(year_label, levels = custom_order))
  
  cld_list <- dunn_results %>%
    group_by(year_label) %>%
    group_split() %>%
    map(function(df) {
      yl      <- as.character(df$year_label[1])
      letters <- get_cld_letters(df)
      if (is.null(letters)) return(NULL)
      tibble(year_label = yl, Zone = names(letters), cld = unname(letters))
    }) %>%
    compact() %>%
    bind_rows() %>%
    mutate(year_label = factor(year_label, levels = custom_order))
  
  # E. ANNOTATIONS
  kw_anno <- kw_results %>%
    filter(p < 0.05) %>%
    mutate(
      year_label = factor(year_label, levels = custom_order),
      Zone  = levels(merged_data$Zone)[1],
      y     = Inf,
      label = if_else(p < 0.001, "KW p < 0.001 ***",
                      paste0("KW p = ", formatC(p, digits = 3), " ", p.signif))
    )
  
  cld_anno <- cld_list %>%
    group_by(year_label) %>%
    filter(n_distinct(cld) > 1) %>%
    ungroup() %>%
    mutate(y = Inf)
  
  # F. PLOT (Boxplots with Zone fill + thin black outline)
  y_axis_lab <- if (grepl("Shannon", plot_label, ignore.case = TRUE)) {
    expression(bold(paste("Shannon Index (", italic("H'"), ")")))
  } else { as.expression(bquote(bold(.(plot_label)))) }
  
  p <- ggplot(merged_data, aes(x = Zone, y = diversity_value)) +
    # Boxplot with Zone background fill, light transparency (0.35), and thin black outline
    geom_boxplot(aes(fill = Zone), color = "black", alpha = 0.35, 
                 outlier.shape = NA, width = 0.6, lwd = 0.3)
  
  if (use_alpha) {
    p <- p + geom_jitter(aes(fill = Zone, alpha = .data[[ALPHA_VAR]]), 
                         shape = 21, color = "white", stroke = 0.35, size = 2.5, width = 0.15)
  } else {
    p <- p + geom_jitter(aes(fill = Zone), 
                         shape = 21, color = "white", stroke = 0.35, size = 2.5, width = 0.15, alpha = 0.85)
  }
  
  p <- p +
    facet_wrap(~year_label, nrow = 1) +
    
    scale_x_discrete(labels = function(x) str_replace_all(x, "_", " ")) +
    
    geom_text(data = kw_anno, aes(x = Zone, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1.3, size = 6, fontface = "italic") +
    
    geom_text(data = cld_anno, aes(x = Zone, y = y, label = cld),
              inherit.aes = FALSE, vjust = 2.9, size = 5, fontface = "bold") +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_npg(name = "Zone") +
    scale_color_npg()
  
  if (use_alpha) {
    alpha_dir    <- if (ALPHA_HIGH_IS_OPAQUE) ALPHA_RANGE else rev(ALPHA_RANGE)
    alpha_vals   <- merged_data[[ALPHA_VAR]]
    alpha_breaks <- pretty(range(alpha_vals, na.rm = TRUE), n = 5)
    
    p <- p + scale_alpha_continuous(
      range  = alpha_dir,
      name   = ALPHA_VAR,
      breaks = alpha_breaks
    )
  }
  
  p <- p +
    labs(x = "Estuary Zone", y = y_axis_lab) +
    
    theme_classic(base_size = 22) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 18),
      axis.title       = element_text(face = "bold", size = 18),
      axis.text.x      = element_text(angle = 45, hjust = 1, color = "black", size = 18),
      axis.text.y      = element_text(color = "black", size = 18),
      plot.title       = element_text(face = "bold", size = 18, hjust = 0),
      legend.position  = "none",
      panel.spacing    = unit(1.5, "lines")
    )
  
  # G. SAVE INDIVIDUAL FILES
  if (save_individual) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    safe_metric <- str_replace_all(metric_pattern, "[^A-Za-z0-9_]", "_")
    
    n_facets   <- n_distinct(merged_data$year_label)
    plot_width <- max(7, n_facets * 2.3)
    
    ggsave(
      filename = file.path(out_dir, paste0("Plot_", marker, "_", safe_metric, ".pdf")),
      plot     = p,
      width    = plot_width,
      height   = 5.5,
      useDingbats = FALSE
    )
    
    ggsave(
      filename = file.path(out_dir, paste0("Plot_", marker, "_", safe_metric, ".png")),
      plot     = p,
      width    = plot_width,
      height   = 5.5
    )
    
    write_csv(kw_results, file.path(out_dir, paste0("Stats_", marker, "_", safe_metric, "_kruskal.csv")))
  }
  
  message("Done.")
  return(invisible(p))
}

# ==============================================================================
# 6. RUN BATCH & COLLECT PLOTS
# ==============================================================================

plot_list <- list()

for (batch in analysis_batches) {
  key <- paste(batch$marker, batch$metric_name, sep = "_")
  result <- try(
    process_and_plot(DATA_DIR, METADATA_PATH, batch$marker, batch$metric_name,
                     batch$label, OUTPUT_DIR, save_individual = TRUE)
  )
  if (!inherits(result, "try-error") && !is.null(result)) {
    plot_list[[key]] <- result
  }
}

# ==============================================================================
# 7. COMBINE INTO ONE PANEL: SHARED X-AXIS, DISTINCT Y-AXES
# ==============================================================================

combine_shared_x <- function(plots) {
  n <- length(plots)
  if (n == 0) return(NULL)
  
  for (i in seq_len(n - 1)) {
    plots[[i]] <- plots[[i]] +
      theme(
        axis.line.x  = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  for (i in setdiff(1:n, c(1, 3))) {
    plots[[i]] <- plots[[i]] +
      theme(strip.text = element_blank())
  }
  
  # wrap_plots(plots, ncol = 1) +
  #   plot_layout(axis_titles = "collect_x")
  
  wrap_plots(plots, ncol = 1) +
    plot_layout(axis_titles = "collect_x") +
    plot_annotation(
      tag_levels = "A" # Automatically assigns A, B, C, D...
    ) & 
    theme(
      plot.tag = element_text(size = 18, face = "bold") # Format the tag size/weight
    )
}

ordered_keys <- c("mifish_shannon", "mifish_faith-pd", "COI_shannon", "COI_faith-pd")
ordered_plots <- plot_list[ordered_keys]
ordered_plots <- compact(ordered_plots)

if (length(ordered_plots) > 0) {
  combined_panel <- combine_shared_x(ordered_plots)
  
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename    = file.path(OUTPUT_DIR, "Combined_Panel_AlphaDiversity.pdf"),
    plot        = combined_panel,
    width       = 13,
    height      = 4.5 * length(ordered_plots),
    useDingbats = FALSE,
    limitsize   = FALSE
  )
  
  ggsave(
    filename  = file.path(OUTPUT_DIR, "Combined_Panel_AlphaDiversity.png"),
    plot      = combined_panel,
    width     = 13,
    height    = 4.5 * length(ordered_plots),
    dpi       = 300,
    limitsize = FALSE
  )
  
  message("\n>>> Combined panel saved to: ", OUTPUT_DIR)
} else {
  warning("No plots were generated — check that input files match the batch patterns.")
}