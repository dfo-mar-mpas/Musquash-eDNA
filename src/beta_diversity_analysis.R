


#!/usr/bin/env Rscript

# ==============================================================================
# BETA DIVERSITY PIPELINE
# PERMANOVA · PERMDISP · envfit · Publication-Quality Visualization
# Alpha transparency driven by environmental variable (e.g. Salinity)
# Combined per-metric figures (patchwork-stacked across markers)
# ==============================================================================

suppressPackageStartupMessages({
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(vegan, ape, ggplot2, dplyr, tidyr, stringr,
                 ggpubr, ggsci, tools, ggrepel, scales, patchwork)
})

# ── NULL-coalescing operator ───────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# CONFIGURATION  — edit only this block
# ==============================================================================

INPUT_DIR   <- "E:/projects/Musquash_data/Analysis/results_diversity_reanalysis/no_prvl_filtering/combined_matrices"
METADATA_IN <- "E:/projects/Musquash_data/Analysis/metadata/comb_station_grp_all_years_metada.tsv"
OUTPUT_DIR  <- "E:/projects/Musquash_data/Analysis/results_diversity_reanalysis/no_prvl_filtering/combined_plots_k8"

FILE_EXT    <- "\\.(tsv|txt|csv)$"
GROUP_COL   <- "Zone"           # Metadata column for grouping / colouring
NMDS_K      <- 6                # NMDS dimensions
SEED        <- 123
PT_SIZE     <- 3.5              # Point size
ELLIPSE_LVL <- 0.95             # Confidence level for stat_ellipse
FONT_SIZE   <- 16               # Base font size (pt) — publication quality

# ── Custom Facet Ordering ─────────────────────────────────────────────────────
CUSTOM_ORDER <- c(
  "2022 NS_MiU", "2023 NS_MiU", "2023 MS_MiM", "2024 MS_MiM", "2025 MS_MiM",
  "2022 NS_LerXT", "2023 NS_LerXT", "2024 MS_Ler", "2025 MS_Ler"
)

# Explicit marker plotting order (uses the RAW marker key as parsed from
# filenames, e.g. "mifish", NOT the display label).
MARKER_ORDER <- c("mifish", "COI")

# Display labels shown on the combined figure — maps raw marker key -> label.
# Everything not listed here falls back to its own raw name.
MARKER_LABELS <- c(mifish = "12S", COI = "COI")

marker_display <- function(marker) {
  lbl <- MARKER_LABELS[marker]
  if (is.na(lbl)) marker else lbl
}

# ── Alpha / transparency settings ─────────────────────────────────────────────
ALPHA_VAR            <- "Salinity" # Metadata column used for transparency
ALPHA_RANGE          <- c(0.30, 0.97)
ALPHA_HIGH_IS_OPAQUE <- TRUE

# ── envfit settings ───────────────────────────────────────────────────────────
ENV_VARS           <- c("Salinity", "Temperature")
ENVFIT_P_CUT       <- 0.05
ENVFIT_SCALE       <- 0.42
ENVFIT_PERMS       <- 999
ENVFIT_COL         <- "grey20"
ENVFIT_LINE_W      <- 0.4
ENVFIT_TXT_SIZE    <- 4.2
SHOW_ENVFIT_ARROWS <- TRUE

# ── Statistical parameters ────────────────────────────────────────────────────
PERM_N <- 999   # Permutations for PERMANOVA and PERMDISP

# ── Plot dimensions ───────────────────────────────────────────────────────────
PLOT_W   <- 14    # inches (single-marker width; combined figure reuses this width)
PLOT_H   <- 5.5   # inches per marker panel
PLOT_DPI <- 600

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# -- Filename parser -----------------------------------------------------------
parse_filename <- function(filepath) {
  fname <- file_path_sans_ext(basename(filepath))
  parts <- str_split(fname, "-")[[1]]
  if (length(parts) >= 4) {
    yr <- parts[2]
    yr_full <- ifelse(nchar(yr) == 2, paste0("20", yr), yr)
    label_val <- paste0(yr_full, " ", parts[3])
    
    data.frame(Filepath = filepath, Filename = fname, Project = parts[1],
               Year = yr_full, FileLabel = label_val, Marker = parts[4], Metric = parts[5],
               stringsAsFactors = FALSE)
  } else NULL
}

# -- Distance matrix loader ----------------------------------------------------
load_dist_matrix <- function(path) {
  d <- tryCatch(
    read.table(path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE),
    error = function(e)
      read.table(path, header = TRUE, row.names = 1, sep = ",", check.names = FALSE)
  )
  as.dist(as.matrix(d))
}

# -- Pairwise PERMANOVA (manual, Bonferroni-corrected) -------------------------
run_pairwise_permanova <- function(D, groups, n_perm = PERM_N, seed = SEED) {
  groups  <- as.character(groups)
  lvls    <- unique(groups)
  if (length(lvls) < 2) return(NULL)
  
  pairs   <- combn(lvls, 2, simplify = TRUE)
  n_pairs <- ncol(pairs)
  
  res <- lapply(seq_len(n_pairs), function(k) {
    g1  <- pairs[1, k]; g2 <- pairs[2, k]
    idx <- which(groups %in% c(g1, g2))
    sub_D   <- as.dist(as.matrix(D)[idx, idx])
    sub_grp <- groups[idx]
    set.seed(seed)
    fit <- adonis2(sub_D ~ sub_grp, permutations = n_perm)
    data.frame(
      Comparison = paste0(g1, " vs ", g2),
      F_value    = round(fit$F[1], 3),
      R2         = round(fit$R2[1], 3),
      p_value    = round(fit$`Pr(>F)`[1], 4),
      stringsAsFactors = FALSE
    )
  })
  
  out        <- do.call(rbind, res)
  out$p_bonf <- round(pmin(out$p_value * n_pairs, 1), 4)
  out
}

# -- Format p-value for labels -------------------------------------------------
fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) return("n.a.")
  if (p < 0.001) "< 0.001" else paste0("= ", round(p, 3))
}

# -- Prepare significant envfit arrows scaled to ordination space --------------
prepare_arrows <- function(ef, ord_df,
                           scale = ENVFIT_SCALE, p_cut = ENVFIT_P_CUT) {
  if (is.null(ef) || is.null(ef$vectors)) return(NULL)
  pv  <- ef$vectors$pvals
  sig <- which(pv <= p_cut)
  if (length(sig) == 0) return(NULL)
  
  sc             <- as.data.frame(scores(ef, "vectors"))[sig, , drop = FALSE]
  sc$Variable    <- rownames(sc)
  sc$p_val       <- pv[sig]
  colnames(sc)[1:2] <- c("Dim1", "Dim2")
  
  max_ord   <- max(abs(c(range(ord_df$Dim1), range(ord_df$Dim2))))
  max_arrow <- max(sqrt(sc$Dim1^2 + sc$Dim2^2))
  if (max_arrow > 0) {
    fac    <- max_ord * scale / max_arrow
    sc$Dim1 <- sc$Dim1 * fac
    sc$Dim2 <- sc$Dim2 * fac
  }
  sc
}

# -- Save a list of data frames to CSV (only if non-empty) --------------------
save_df <- function(lst, filepath) {
  df <- do.call(rbind, Filter(Negate(is.null), lst))
  if (!is.null(df) && nrow(df) > 0)
    write.csv(df, filepath, row.names = FALSE)
}

# ==============================================================================
# SETUP
# ==============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

meta <- read.table(METADATA_IN, header = TRUE, sep = "\t",
                   check.names = FALSE, comment.char = "")

sid_idx <- which(tolower(colnames(meta)) %in%
                   c("sampleid", "sample-id", "sample_id", "id"))[1]
colnames(meta)[sid_idx] <- "SampleID"

avail_env <- ENV_VARS[ENV_VARS %in% colnames(meta)]
if (length(avail_env) < length(ENV_VARS)) {
  warning("envfit vars not found in metadata (will skip): ",
          paste(setdiff(ENV_VARS, avail_env), collapse = ", "))
}

files        <- list.files(INPUT_DIR, pattern = FILE_EXT, full.names = TRUE)
design_table <- do.call(rbind, lapply(files, parse_filename))
combos       <- design_table %>% distinct(Marker, Metric)

# -- Master accumulators for all markers, metrics, and years --
master_perm_glob <- list()
master_perm_pair <- list()
master_disp_glob <- list()
master_disp_pair <- list()
master_pcoa_env  <- list()
master_nmds_env  <- list()

# -- NEW: accumulators for combined per-metric figures --
# Each is a named list keyed by Metric; each element is itself a named list
# keyed by Marker holding the finished ggplot object for that marker.
pcoa_plots_by_metric <- list()
nmds_plots_by_metric <- list()

# ==============================================================================
# MAIN LOOP  — iterate over each Marker × Metric combination
# ==============================================================================

for (i in seq_len(nrow(combos))) {
  
  curr_marker <- combos$Marker[i]
  curr_metric <- combos$Metric[i]
  tag         <- paste0(curr_marker, "_", curr_metric)
  
  cat(sprintf("\n══════════════════════════════════════════\n"))
  cat(sprintf("  Marker: %s  |  Metric: %s\n", curr_marker, curr_metric))
  cat(sprintf("══════════════════════════════════════════\n"))
  
  year_files <- design_table %>%
    filter(Marker == curr_marker, Metric == curr_metric) %>%
    arrange(Year)
  
  # -- Per-year accumulators --------------------------------------------------
  pcoa_pts  <- list(); pcoa_env  <- list()
  nmds_pts  <- list(); nmds_env  <- list()
  
  perm_glob_all <- list(); perm_pair_all <- list()
  disp_glob_all <- list(); disp_pair_all <- list()
  pcoa_env_stats <- list()
  nmds_env_stats <- list()
  
  # ── YEAR LOOP ─────────────────────────────────────────────────────────────
  for (j in seq_len(nrow(year_files))) {
    
    f_path <- year_files$Filepath[j]
    f_year <- year_files$Year[j]
    cat(sprintf("\n  ▸ Year: %s\n", f_year))
    
    D            <- load_dist_matrix(f_path)
    used_samples <- labels(D)
    
    sub_meta <- meta %>%
      filter(SampleID %in% used_samples) %>%
      arrange(match(SampleID, used_samples))
    
    if (nrow(sub_meta) < 3) {
      cat("    Skipping – fewer than 3 matched samples.\n"); next
    }
    
    groups <- as.factor(sub_meta[[GROUP_COL]])
    
    env_df <- sub_meta %>%
      select(all_of(avail_env)) %>%
      mutate(across(everything(), as.numeric))
    
    # ─────────────────────────────────────────────────────────────────────
    # 1. PERMANOVA
    # ─────────────────────────────────────────────────────────────────────
    set.seed(SEED)
    perm_fit <- adonis2(D ~ groups, permutations = PERM_N)
    perm_p   <- perm_fit$`Pr(>F)`[1]
    perm_F   <- round(perm_fit$F[1], 3)
    perm_R2  <- round(perm_fit$R2[1], 3)
    
    master_perm_glob[[length(master_perm_glob) + 1]] <- data.frame(
      Marker = curr_marker, Metric = curr_metric, Year = f_year,
      F_value = perm_F, R2 = perm_R2, p_value = round(perm_p, 4)
    )
    
    pw_perm <- run_pairwise_permanova(D, groups)
    if (!is.null(pw_perm)) {
      pw_perm$Marker <- curr_marker
      pw_perm$Metric <- curr_metric
      pw_perm$Year   <- f_year
      master_perm_pair[[length(master_perm_pair) + 1]] <- pw_perm
      cat("    Pairwise PERMANOVA:\n")
      print(pw_perm, row.names = FALSE)
    }
    
    # ─────────────────────────────────────────────────────────────────────
    # 2. PERMDISP
    # ─────────────────────────────────────────────────────────────────────
    disp_obj  <- betadisper(D, groups, type = "median")
    set.seed(SEED)
    disp_test <- permutest(disp_obj, permutations = PERM_N)
    disp_p    <- disp_test$tab$`Pr(>F)`[1]
    disp_F    <- round(disp_test$tab$F[1], 3)
    
    master_disp_glob[[length(master_disp_glob) + 1]] <- data.frame(
      Marker = curr_marker, Metric = curr_metric, Year = f_year,
      F_value = disp_F, p_value = round(disp_p, 4)
    )
    
    if (!is.na(disp_p) && disp_p < 0.05 && length(levels(groups)) > 2) {
      tuk            <- as.data.frame(TukeyHSD(disp_obj)$group)
      tuk$Comparison <- rownames(tuk)
      tuk$Marker     <- curr_marker
      tuk$Metric     <- curr_metric
      tuk$Year       <- f_year
      master_disp_pair[[length(master_disp_pair) + 1]] <- tuk
      cat("    Pairwise PERMDISP (TukeyHSD):\n")
      print(tuk[, c("Comparison", "diff", "lwr", "upr", "p adj")], row.names = FALSE)
    }
    
    # ─────────────────────────────────────────────────────────────────────
    # 3. Build facet strip label  (year + ordination info + test p-values)
    # ─────────────────────────────────────────────────────────────────────
    f_label <- year_files$FileLabel[j]
    
    make_strip <- function(extra_line) {
      paste0(f_label, "\n", extra_line)
    }
    
    # ─────────────────────────────────────────────────────────────────────
    # 4. PCoA
    # ─────────────────────────────────────────────────────────────────────
    pcoa_res <- pcoa(D, correction = "cailliez")
    
    pc_raw          <- as.data.frame(pcoa_res$vectors[, 1:2])
    colnames(pc_raw) <- c("Dim1", "Dim2")
    pc_raw$Dim1     <- pc_raw$Dim1 - mean(pc_raw$Dim1)
    pc_raw$Dim2     <- pc_raw$Dim2 - mean(pc_raw$Dim2)
    pc_raw$SampleID <- rownames(pcoa_res$vectors)
    
    eig     <- pcoa_res$values$Rel_corr_eig %||% pcoa_res$values$Eigenvalues
    var_exp <- round(eig / sum(eig) * 100, 1)
    
    fl_pcoa <- make_strip(sprintf("PC1: %.1f%%  |  PC2: %.1f%%",
                                  var_exp[1], var_exp[2]))
    
    pc_df              <- left_join(pc_raw, sub_meta, by = "SampleID")
    pc_df$FacetLabel   <- fl_pcoa
    pc_df$RealYear     <- f_year
    pcoa_pts[[j]]      <- pc_df
    
    if (nrow(env_df) == nrow(pc_df) && length(avail_env) > 0) {
      ef_pc <- tryCatch(
        envfit(as.matrix(pc_df[, c("Dim1", "Dim2")]), env_df,
               permutations = ENVFIT_PERMS, na.rm = TRUE),
        error = function(e) { cat("    envfit (PCoA) failed:", e$message, "\n"); NULL }
      )
      
      if (!is.null(ef_pc) && !is.null(ef_pc$vectors)) {
        master_pcoa_env[[length(master_pcoa_env) + 1]] <- data.frame(
          Marker = curr_marker, Metric = curr_metric, Year = f_year,
          Variable = rownames(ef_pc$vectors$arrows),
          R2 = ef_pc$vectors$r,
          p_value = ef_pc$vectors$pvals,
          stringsAsFactors = FALSE
        )
      }
      
      arr_pc <- prepare_arrows(ef_pc, pc_df)
      if (!is.null(arr_pc)) {
        arr_pc$FacetLabel <- fl_pcoa
        pcoa_env[[j]]     <- arr_pc
        cat(sprintf("    envfit (PCoA) – %d significant vector(s).\n", nrow(arr_pc)))
      }
    }
    
    # ─────────────────────────────────────────────────────────────────────
    # 5. NMDS
    # ─────────────────────────────────────────────────────────────────────
    set.seed(SEED)
    nmds_res <- try(
      metaMDS(D, k = NMDS_K, trymax = 1000,
              maxit = 1000,trace = 0,
              autotransform = FALSE, wascores = FALSE),
      silent = TRUE
    )
    
    if (!inherits(nmds_res, "try-error")) {
      
      if (!nmds_res$converged) {
        cat(sprintf("    [NOTE] NMDS did not reach strict convergence. Using best solution found (stress = %.4f).\n",
                    nmds_res$stress))
        fl_nmds <- make_strip(sprintf("Stress: %.3f (nc)", nmds_res$stress))
      } else {
        cat(sprintf("    NMDS        : stress = %.4f\n", nmds_res$stress))
        fl_nmds <- make_strip(sprintf("Stress: %.3f", nmds_res$stress))
      }
      
      # nm_raw           <- as.data.frame(scores(nmds_res, display = "sites"))
      # colnames(nm_raw) <- c("Dim1", "Dim2")
      # nm_raw$SampleID  <- rownames(nm_raw)
      # 
      # nm_df            <- left_join(nm_raw, sub_meta, by = "SampleID")
      nm_scores <- scores(nmds_res, display = "sites")
      
      # Explicitly take only the first two NMDS axes for plotting
      nm_raw <- data.frame(
        Dim1 = nm_scores[, 1],
        Dim2 = nm_scores[, 2],
        SampleID = rownames(nm_scores),
        check.names = FALSE
      )
      
      nm_df <- left_join(nm_raw, sub_meta, by = "SampleID")
      nm_df$FacetLabel <- fl_nmds
      nm_df$RealYear   <- f_year
      nmds_pts[[j]]    <- nm_df
      
      if (nrow(env_df) == nrow(nm_df) && length(avail_env) > 0) {
        ef_nm <- tryCatch(
          envfit(nmds_res, env_df, permutations = ENVFIT_PERMS, na.rm = TRUE),
          error = function(e) { cat("    envfit (NMDS) failed:", e$message, "\n"); NULL }
        )
        
        if (!is.null(ef_nm) && !is.null(ef_nm$vectors)) {
          master_nmds_env[[length(master_nmds_env) + 1]] <- data.frame(
            Marker = curr_marker, Metric = curr_metric, Year = f_year,
            Variable = rownames(ef_nm$vectors$arrows),
            R2 = ef_nm$vectors$r,
            p_value = ef_nm$vectors$pvals,
            stringsAsFactors = FALSE
          )
        }
        
        arr_nm <- prepare_arrows(ef_nm, nm_df)
        if (!is.null(arr_nm)) {
          arr_nm$FacetLabel <- fl_nmds
          nmds_env[[j]]     <- arr_nm
          cat(sprintf("    envfit (NMDS) – %d significant vector(s).\n", nrow(arr_nm)))
        }
      }
      
    } else {
      cat("    NMDS threw a fatal error — skipping.\n")
    }
    
  }  # ── end year loop ──────────────────────────────────────────────────────
  
  # ── PLOT GENERATION FUNCTION ─────────────────────────────────────────────
  # UNCHANGED — same facet strips, envfit arrows, alpha-by-salinity,
  # custom facet ordering, legends, and theme as before.
  make_plot <- function(pts_list, env_list, method_name) {
    
    pts_list <- Filter(Negate(is.null), pts_list)
    if (length(pts_list) == 0) return(NULL)
    
    big <- do.call(rbind, pts_list)
    
    facet_levels <- sapply(CUSTOM_ORDER, function(lbl) {
      match_idx <- which(grepl(paste0("^", lbl, "(\\n|$)"), unique(big$FacetLabel)))
      if (length(match_idx) > 0) unique(big$FacetLabel)[match_idx] else NULL
    })
    facet_levels <- unlist(Filter(Negate(is.null), facet_levels))
    
    big$FacetLabel <- factor(big$FacetLabel, levels = facet_levels)
    
    use_alpha <- ALPHA_VAR %in% colnames(big) &&
      any(!is.na(as.numeric(big[[ALPHA_VAR]])))
    if (use_alpha) big[[ALPHA_VAR]] <- as.numeric(big[[ALPHA_VAR]])
    
    arrows_df  <- NULL
    env_valid  <- Filter(Negate(is.null), env_list)
    if (length(env_valid) > 0) {
      arrows_df             <- do.call(rbind, env_valid)
      arrows_df$FacetLabel  <- factor(arrows_df$FacetLabel, levels = levels(big$FacetLabel))
      arrows_df$Variable    <- substr(arrows_df$Variable, 1, 1)
    }
    
    p <- ggplot(big, aes(x = Dim1, y = Dim2))
    
    p <- p + stat_ellipse(
      aes(color = .data[[GROUP_COL]]),
      level = ELLIPSE_LVL, linetype = "dashed",
      linewidth = 0.4, show.legend = FALSE
    )
    
    if (use_alpha) {
      p <- p + geom_point(
        aes(fill  = .data[[GROUP_COL]],
            alpha = .data[[ALPHA_VAR]]),
        shape = 21, color = "white", stroke = 0.35, size = PT_SIZE
      )
    } else {
      p <- p + geom_point(
        aes(fill = .data[[GROUP_COL]]),
        shape = 21, color = "white", stroke = 0.35,
        size = PT_SIZE, alpha = 0.82
      )
    }
    
    if (!is.null(arrows_df) && SHOW_ENVFIT_ARROWS) {
      p <- p +
        geom_segment(
          data = arrows_df,
          aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
          arrow = arrow(length = unit(0.14, "cm"), type = "open", angle = 22),
          colour = ENVFIT_COL,
          linewidth = ENVFIT_LINE_W,
          inherit.aes = FALSE
        ) +
        geom_text_repel(
          data        = arrows_df,
          aes(x = Dim1, y = Dim2, label = Variable),
          colour      = "grey35",
          size        = ENVFIT_TXT_SIZE,
          box.padding = 0.35,
          max.overlaps = 25,
          inherit.aes = FALSE
        )
    }
    
    p <- p +
      facet_wrap(~FacetLabel, scales = "free", nrow = 1) +
      scale_fill_npg(name = GROUP_COL) +
      scale_color_npg()
    
    if (use_alpha) {
      alpha_dir    <- if (ALPHA_HIGH_IS_OPAQUE) ALPHA_RANGE else rev(ALPHA_RANGE)
      alpha_vals   <- big[[ALPHA_VAR]]
      alpha_breaks <- pretty(range(alpha_vals, na.rm = TRUE), n = 5)
      
      p <- p + scale_alpha_continuous(
        range  = alpha_dir,
        name   = ALPHA_VAR,
        breaks = alpha_breaks,
        guide  = guide_legend(
          title.position = "top",
          title.hjust    = 0.5,
          nrow           = length(alpha_breaks),
          order          = 2,
          override.aes   = list(
            fill   = "grey42",
            color  = "white",
            shape  = 21,
            size   = 3.8,
            stroke = 0.3
          )
        )
      )
    }
    
    p <- p + guides(
      fill = guide_legend(
        title.position = "top",
        title.hjust    = 0.5,
        order          = 1,
        override.aes   = list(size = 4, stroke = 0.3, color = "white", alpha = 0.92)
      )
    )
    
    p <- p +
      labs(
        x = paste(method_name, "Axis 1"),
        y = paste(method_name, "Axis 2")
      ) +
      theme_classic(base_size = FONT_SIZE) +
      theme(
        # strip.background = element_rect(fill = "grey94", colour = NA),
        strip.background = element_blank(),
        strip.text       = element_text(
          face = "bold", size = FONT_SIZE - 3,
          lineheight = 1.2,
          margin = margin(t = 6, b = 6, l = 4, r = 4)
        ),
        axis.title  = element_text(face = "bold", size = FONT_SIZE),
        axis.text   = element_text(size = FONT_SIZE - 2, colour = "grey15"),
        axis.line   = element_line(linewidth = 0.55, colour = "grey25"),
        axis.ticks  = element_line(linewidth = 0.45, colour = "grey25"),
        legend.title     = element_text(face = "bold", size = FONT_SIZE - 2),
        legend.text      = element_text(size = FONT_SIZE - 3),
        legend.key.size  = unit(0.52, "cm"),
        legend.spacing.y = unit(0.15, "cm"),
        legend.position  = "right",
        legend.box       = "vertical",
        legend.box.spacing = unit(0.3, "cm"),
        panel.spacing  = unit(1.2, "lines"),
        plot.margin = margin(t = 8, r = 14, b = 8, l = 8)
      )
    
    p
  }
  
  # ── BUILD PLOTS, BUT STASH INSTEAD OF SAVING IMMEDIATELY ────────────────
  p_pcoa <- make_plot(pcoa_pts, pcoa_env, "PCoA")
  if (!is.null(p_pcoa)) {
    if (is.null(pcoa_plots_by_metric[[curr_metric]])) {
      pcoa_plots_by_metric[[curr_metric]] <- list()
    }
    pcoa_plots_by_metric[[curr_metric]][[curr_marker]] <- p_pcoa
    cat(sprintf("\n  ✓ PCoA plot generated (stashed): %s / %s\n", curr_marker, curr_metric))
  }
  
  p_nmds <- make_plot(nmds_pts, nmds_env, "NMDS")
  if (!is.null(p_nmds)) {
    if (is.null(nmds_plots_by_metric[[curr_metric]])) {
      nmds_plots_by_metric[[curr_metric]] <- list()
    }
    nmds_plots_by_metric[[curr_metric]][[curr_marker]] <- p_nmds
    cat(sprintf("  ✓ NMDS plot generated (stashed): %s / %s\n", curr_marker, curr_metric))
  }
  
  cat(sprintf("  ✓ Stats accumulated for: %s\n", tag))
  
}  # ── end marker × metric loop ─────────────────────────────────────────────


# ==============================================================================
# COMBINE PER-METRIC FIGURES WITH PATCHWORK  (one PDF per metric, per method)
# ==============================================================================
# Each marker's plot already carries its own facet strips, legend, envfit
# arrows, and alpha-by-salinity encoding — those are untouched. Here we only
# stack the marker panels for a given metric into a single tall figure, and
# collect legends so they appear once.

combine_and_save <- function(plots_by_metric, method_name, out_prefix) {
  for (metric in names(plots_by_metric)) {
    marker_plots <- plots_by_metric[[metric]]
    marker_plots <- Filter(Negate(is.null), marker_plots)
    if (length(marker_plots) == 0) next
    
    # Stack vertically in marker order (order of first appearance in combos)
    # marker_order <- combos %>% filter(Metric == metric) %>% pull(Marker) %>% unique()
    # marker_order <- marker_order[marker_order %in% names(marker_plots)]
    # marker_plots <- marker_plots[marker_order]
    marker_order <- MARKER_ORDER[MARKER_ORDER %in% names(marker_plots)]
    # Fallback: append any markers not covered by MARKER_ORDER, preserving
    # their natural (combos) order, so nothing silently gets dropped.
    leftover <- setdiff(names(marker_plots), marker_order)
    if (length(leftover) > 0) {
      leftover <- combos %>% filter(Metric == metric, Marker %in% leftover) %>%
        pull(Marker) %>% unique()
      marker_order <- c(marker_order, leftover)
    }
    marker_plots <- marker_plots[marker_order]
    
    combined <- wrap_plots(marker_plots, ncol = 1) +
      plot_layout(guides = "collect") +
      # plot_annotation(tag_levels = "A")
      plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag = element_text(face = "bold")
      )
      )
    
    n_markers <- length(marker_plots)
    out_path  <- file.path(OUTPUT_DIR, paste0(out_prefix, "_", metric, "_combined.pdf"))
    
    ggsave(
      out_path, combined,
      width = PLOT_W, height = PLOT_H * n_markers, dpi = PLOT_DPI,
      limitsize = FALSE
    )
    cat(sprintf("  ✓ Combined %s plot saved for metric '%s': %s\n",
                method_name, metric, basename(out_path)))
  }
}

combine_and_save(pcoa_plots_by_metric, "PCoA", "PCoA")
combine_and_save(nmds_plots_by_metric, "NMDS", "NMDS")


# ── SAVE MASTER STATISTICAL RESULTS ─────────────────────────────────────────
save_df(master_perm_glob, file.path(OUTPUT_DIR, "PERMANOVA_Global_All.csv"))
save_df(master_perm_pair, file.path(OUTPUT_DIR, "PERMANOVA_Pairwise_All.csv"))
save_df(master_disp_glob, file.path(OUTPUT_DIR, "PERMDISP_Global_All.csv"))
save_df(master_disp_pair, file.path(OUTPUT_DIR, "PERMDISP_Pairwise_All.csv"))
save_df(master_pcoa_env,  file.path(OUTPUT_DIR, "envfit_stats_PCoA_All.csv"))
save_df(master_nmds_env,  file.path(OUTPUT_DIR, "envfit_stats_NMDS_All.csv"))

cat("\n══════════════════════════════════════════\n")
cat("  Pipeline complete!\n")
cat("══════════════════════════════════════════\n")
