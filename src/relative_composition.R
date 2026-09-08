# ============================================================
#  Publication-Grade Community Composition Visualizations
#  NPG Style | Multi-Year, Multi-Zone Taxonomic Data
#  Target: Q1 Journal (Nature / Cell / Science family)
# ============================================================

# ── 0. PACKAGES ─────────────────────────────────────────────
required_pkgs <- c(
  "tidyverse", "ggsci", "ggalluvial", "scales",
  "ggtext", "patchwork", "RColorBrewer", "forcats", "stringr", "rfishbase"
)
new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ── 1. USER SETTINGS ────────────────────────────────────────
INPUT_FILE   <- "E:/projects/Musquash_data/Analysis/results2/mifish_species_viz_filtered/Musq-mifish_zones_years.csv"
OUTPUT_DIR   <- "E:/projects/Musquash_data/Analysis/results2/mifish_species_viz_filtered/Musq-mifish-relative-composition2"
# INPUT_FILE   <- "E:/projects/Musquash_data/Analysis/results2/COI_species_viz_filtered/Musq-COI_zones_years.csv"
# OUTPUT_DIR   <- "E:/projects/Musquash_data/Analysis/results2/COI_species_viz_filtered/Musq-COI-relative-composition_23incl"
TOP_N        <- 20
MIN_ABUND    <- 0
DPI          <- 350
BASE_FONT    <- 15

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ── 2. LOAD DATA ─────────────────────────────────────────────
df_raw <- read.csv(INPUT_FILE, row.names = 1, check.names = FALSE)
stopifnot(nrow(df_raw) > 0, ncol(df_raw) > 0)

# ── 3. PARSE ROW NAMES → TAXONOMY TABLE ─────────────────────
parse_lineage <- function(lin) {
  extract <- function(tag) {
    str_extract(lin, paste0("(?<=", tag, "__)[^;]+")) %>%
      str_trim() %>% str_replace_all("_", " ") %>%
      na_if("") %>% { if (is.na(.)) "Unknown" else . }
  }
  kingdom <- extract("k")
  phylum  <- extract("p")
  class_  <- extract("c")
  order_  <- extract("o")
  family_ <- extract("f")
  genus_  <- extract("g")
  
  sp_block <- str_extract(lin, "s__[^;]*$") %>%
    str_replace("^s__", "") %>% str_trim()
  
  common <- str_extract(sp_block, "(?<=\\()[^)]+(?=\\))")
  sci    <- str_replace(sp_block, "\\s*\\([^)]*\\)", "") %>%
    str_trim() %>% str_replace_all("_", " ")
  
  sci_label    <- ifelse(is.na(sci)    | sci    == "", genus_, sci)
  common_label <- ifelse(is.na(common) | common == "", sci_label, common)
  
  tibble(
    Lineage      = lin,
    Kingdom      = kingdom,
    Phylum       = phylum,
    Class        = class_,
    Order        = order_,
    Family       = family_,
    Genus        = genus_,
    Species      = sci_label,
    CommonName   = common_label,
    # DisplayLabel = ifelse(
    #   !is.na(common) & common != "" & common != sci_label,
    #   paste0(sci_label, "\n(", common, ")"),
    #   sci_label
    # )
    DisplayLabel = ifelse(
      !is.na(common) & common != "" & common != sci_label,
      paste0("<i>", sci_label, "</i><br><b>(", common, ")</b>"),
      ifelse(sci_label == "Other", "Other", paste0("<i>", sci_label, "</i>"))
    )
  )
}

tax_tbl <- map_dfr(rownames(df_raw), parse_lineage)

# # ── 3.5 RETRIEVE ECOLOGICAL TRAITS FROM FISHBASE ──────────────
# cat("Retrieving ecological traits from FishBase...\n")
# 
# unique_spp <- unique(tax_tbl$Species)
# fb_spp     <- rfishbase::species(unique_spp)
# 
# trait_tbl <- fb_spp %>%
#   select(Species, Fresh, Brack, Saltwater, DemersPelag) %>%
#   mutate(across(c(Fresh, Brack, Saltwater), ~ if_else(!is.na(.) & . != 0, 1, 0))) %>%
#   mutate(
#     salinity_class = case_when(
#       Fresh == 1 & Brack == 0 & Saltwater == 0 ~ "Freshwater",
#       Fresh == 0 & Brack == 1 & Saltwater == 0 ~ "Brackish",
#       Fresh == 0 & Brack == 0 & Saltwater == 1 ~ "Marine",
#       Brack == 1 & (Fresh == 1 | Saltwater == 1) ~ "Euryhaline",
#       TRUE ~ "Unknown"
#     ),
#     Habitat = if_else(is.na(DemersPelag), "Unknown", DemersPelag)
#   ) %>%
#   select(Species, salinity_class, Habitat) %>%
#   distinct(Species, .keep_all = TRUE)
# 
# tax_tbl <- tax_tbl %>%
#   left_join(trait_tbl, by = "Species") %>%
#   mutate(
#     salinity_class = replace_na(salinity_class, "Unknown"),
#     Habitat        = replace_na(Habitat, "Unknown")
#   )

# ── 4. PARSE COLUMN NAMES → SAMPLE METADATA ─────────────────
parse_colname <- function(cn) {
  yr   <- str_extract(cn, "\\b(19|20)\\d{2}\\b")
  zone <- str_replace(cn, paste0("[\\-_]?", yr, ".*"), "") %>% str_trim()
  lbl  <- str_replace(cn, paste0(".*", yr, "[\\-_]?"), "") %>% str_trim()
  lbl  <- if (is.na(lbl) | lbl == "") "—" else lbl
  tibble(
    ColName   = cn,
    Zone      = zone,
    Year      = as.integer(yr),
    Label     = lbl,
    YearLabel = paste0(yr, if (lbl != "—") paste0("\n", lbl) else "")
  )
}

meta_tbl <- map_dfr(colnames(df_raw), parse_colname)

yr_order <- meta_tbl %>% arrange(Year, Label) %>% pull(YearLabel) %>% unique()
zo_order <- sort(unique(meta_tbl$Zone))

meta_tbl <- meta_tbl %>%
  mutate(
    YearLabel = factor(YearLabel, levels = yr_order),
    Zone      = factor(Zone,      levels = zo_order)
  )

# ── 5. RESHAPE TO LONG FORMAT & COMPUTE RELATIVE ABUNDANCE ──
df_long <- df_raw %>%
  rownames_to_column("Lineage") %>%
  pivot_longer(-Lineage, names_to = "ColName", values_to = "ReadCount") %>%
  left_join(tax_tbl,  by = "Lineage") %>%
  left_join(meta_tbl, by = "ColName") %>%
  group_by(ColName) %>%
  mutate(
    TotalReads = sum(ReadCount, na.rm = TRUE),
    RelAbund   = if_else(TotalReads > 0, ReadCount / TotalReads * 100, 0)
  ) %>%
  ungroup()

# ── 6. BUILD GLOBAL COLOUR PALETTE ──────────────────────────
# npg10   <- pal_npg("nrc")(10)
npg_cols    <- pal_npg("nrc")(10)
jama_cols   <- pal_jama("default")(7)
nejm_cols   <- pal_nejm("default")(8)
lancet_cols <- pal_lancet("lanonc")(9)

# extra14 <- c(
#   "#7E6148", "#B09C85", "#00A087", "#3C5488",
#   "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
#   "#7CAE00", "#C77CFF", "#00BFC4", "#F8766D",
#   "#619CFF", "#00BA38"
# )
# palette_pool <- unique(c(npg10, extra14))
journal_pool <- unique(c(npg_cols, jama_cols, nejm_cols, lancet_cols))

# 3. Bring in your custom choices, making sure no duplicates remain
extra_custom <- c("#7E6148", "#B09C85", "#7CAE00", "#C77CFF", "#00BFC4", "#F8766D", "#619CFF", "#00BA38")

palette_pool <- unique(c(journal_pool, extra_custom))

species_rank <- df_long %>%
  group_by(Species, DisplayLabel) %>%
  summarise(MeanAbund = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MeanAbund))

top_species  <- head(species_rank$Species, TOP_N)
label_lookup <- setNames(species_rank$DisplayLabel, species_rank$Species)

sp_colors <- setNames(palette_pool[seq_along(top_species)], top_species)
sp_colors["Other"] <- "#BBBBBB"

df_long <- df_long %>%
  mutate(
    SpeciesGrp   = if_else(Species %in% top_species, Species, "Other"),
    SpeciesGrp   = factor(SpeciesGrp, levels = c(top_species, "Other")),
    DisplayLabel = if_else(Species %in% top_species,
                           label_lookup[Species], "Other"),
    DisplayLabel = factor(DisplayLabel,
                          levels = c(label_lookup[top_species], "Other"))
  )

# ── 7. SHARED NPG THEME ──────────────────────────────────────
theme_npg_pub <- function(base_size = BASE_FONT) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text              = element_text(family = "sans", colour = "#1A1A1A"),
      plot.title        = element_text(size = base_size + 4, face = "bold",
                                       hjust = 0, margin = margin(b = 4)),
      plot.subtitle     = element_text(size = base_size + 1, hjust = 0,
                                       colour = "#555555", margin = margin(b = 10)),
      plot.caption      = element_text(size = base_size - 2, hjust = 1,
                                       colour = "#888888", margin = margin(t = 8)),
      axis.title        = element_text(size = base_size, face = "bold"),
      axis.text         = element_text(size = base_size - 2, colour = "#1A1A1A"),
      axis.line         = element_line(colour = "#333333", linewidth = 0.45),
      axis.ticks        = element_line(colour = "#333333", linewidth = 0.35),
      axis.ticks.length = unit(2.5, "pt"),
      strip.background  = element_rect(fill = "#2B547E", colour = NA),
      strip.text        = element_text(colour = "white", face = "bold",
                                       size = base_size - 0.5, margin = margin(3, 4, 3, 4)),
      strip.clip        = "off",
      panel.grid.major.y = element_line(colour = "#EBEBEB", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(0.45, "cm"),
      legend.position    = "right",
      
      # BUMPED UP FOR Q1 VISIBILITY:
      legend.title       = element_text(face = "bold", size = base_size + 2), # Increased from -1
      # legend.text        = element_text(size = base_size + 1),               # Increased from -2
      # legend.text        = ggtext::element_markdown(size = base_size + 1, hjust = 0),
      legend.text        = ggtext::element_markdown(
        size = base_size + 1, 
        hjust = 0,
        margin = margin(t = 1, b = 1, l = 2, unit = "mm") 
      ),
      legend.key.size    = unit(0.75, "cm"),                                 # Made keys slightly larger to match text

      legend.key         = element_rect(fill = NA, colour = NA),
      legend.frame       = element_blank(),
      legend.background  = element_blank(),
      legend.spacing.y = unit(6, "mm"),
      legend.spacing.x = unit(6, "mm"),
      legend.key.height = unit(14, "mm"),
      legend.key.width  = unit(6, "mm"),
      legend.box.spacing = unit(6, "mm"),
      plot.margin        = margin(12, 12, 10, 10),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.background   = element_rect(fill = "white", colour = NA)
      
      
      # 
      # legend.text        = ggtext::element_markdown(
      #   size = base_size + 1, 
      #   hjust = 0,
      #   margin = margin(t = 6, b = 6, l = 4, unit = "mm") 
      # ),
      # 
      # legend.key.size    = unit(1.5, "cm"),                                 
      # legend.key         = element_rect(fill = NA, colour = NA),
      # legend.frame       = element_blank(),
      # legend.background  = element_blank(),
      # 
      # # INCREASE SPACING.Y (Must be used in combination with byrow = TRUE in your guides)
      # legend.spacing.y   = unit(16, "mm"), 
      # legend.spacing.x   = unit(6, "mm"),
    )
}



scale_fill_species <- function() {
  scale_fill_manual(
    values = sp_colors,
    name   = "Species",
    limits = c(top_species, "Other"),   
    drop   = FALSE,                     
    labels = c(label_lookup[top_species], "Other"),
    guide  = guide_legend(
      ncol         = 1,
      byrow        = TRUE,               # 1. ADD THIS LINE
      keyheight    = unit(1, "cm"),   # 2. MATCHED TO YOUR NEW KEY SIZE
      keywidth     = unit(0.55, "cm"),   
      override.aes = list(alpha = 1)
    )
  )
}

scale_color_species <- function() {
  scale_color_manual(
    values = sp_colors,
    name   = "Species",
    labels = c(label_lookup[top_species], "Other"),
    guide  = guide_legend(
      ncol         = 1,
      byrow        = TRUE,               # 3. ADD THIS LINE
      keyheight    = unit(1, "cm"),   # 4. MATCHED TO YOUR NEW KEY SIZE
      keywidth     = unit(0.55, "cm"),   
      override.aes = list(size = 4, alpha = 1)
    )
  )
}

# # ============================================================
# #  PLOT 1 — FACETED STACKED BAR PLOT
# # ============================================================
# cat("Building Plot 1: Faceted Stacked Bar ...\n")
# 
# # FIX: sum within each sample first, then average across samples.
# # The old code took mean(RelAbund) directly across all species rows
# # sharing a SpeciesGrp, which under-counted multi-species groups
# # (especially "Other") and caused bars to fall well below 100%.
# df_bar <- df_long %>%
#   group_by(ColName, Zone, YearLabel, SpeciesGrp, DisplayLabel) %>%
#   summarise(RA = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%   # ← sum per sample
#   group_by(Zone, YearLabel, SpeciesGrp, DisplayLabel) %>%
#   summarise(RA = mean(RA, na.rm = TRUE), .groups = "drop")            # ← then average
# 
# p1 <- ggplot(df_bar,
#              aes(x = Zone, y = RA, fill = SpeciesGrp)) +
#   geom_col(position = "stack", width = 0.72,
#            colour = "white", linewidth = 0.18) +
#   facet_grid(
#     cols = vars(YearLabel),
#     switch = "x"
#   ) +
#   scale_fill_species() +
#   scale_y_continuous(
#     labels   = function(x) paste0(x, "%"),
#     expand   = expansion(mult = c(0, 0.02)),
#     limits   = c(0, 100)
#   ) +
#   labs(
#     title    = "Spatio Temporal Trends - Community Composition",
#     subtitle = "Mean relative abundance (%) by zone and sampling period",
#     x        = NULL,
#     y        = "Relative Abundance (%)",
#     caption  = paste0("Top ", TOP_N, " taxa by mean abundance shown individually; remainder pooled as 'Other'.")
#   ) +
#   theme_npg_pub() +
#   theme(
#     axis.text.x    = element_text(angle = 35, hjust = 1, size = BASE_FONT - 1.5),
#     strip.placement = "outside"
#   )
# 
# out1_pdf <- file.path(OUTPUT_DIR, "Fig1_stacked_bar.pdf")
# out1_png <- file.path(OUTPUT_DIR, "Fig1_stacked_bar.png")
# ggsave(out1_pdf, p1,
#        width  = max(12, length(yr_order) * 2.2 + 3),
#        height = max(6,  length(zo_order) * 1.5 + 2),
#        device = cairo_pdf)
# ggsave(out1_png, p1,
#        width  = max(12, length(yr_order) * 2.2 + 3),
#        height = max(6,  length(zo_order) * 1.5 + 2),
#        dpi = DPI)
# cat("  ✓ Saved:", out1_pdf, "&", out1_png, "\n")
# 
# # ============================================================
# #  PLOT 2 — FACETED BUBBLE / DOT PLOT
# # ============================================================
# cat("Building Plot 2: Faceted Bubble Plot ...\n")
# 
# df_bubble <- df_long %>%
#   filter(RelAbund > 0) %>%
#   group_by(Zone, YearLabel, Species, DisplayLabel, salinity_class, Habitat, SpeciesGrp) %>%
#   summarise(RA = mean(RelAbund, na.rm = TRUE), .groups = "drop")
# 
# salinity_colors <- c(
#   "Freshwater" = "#00A087",
#   "Euryhaline"   = "#4DBBD5",
#   "Marine"     = "#3C5488",
#   "Brackish" = "#91D1C2",
#   "Unknown"    = "#B09C85"
# )
# 
# sp_ordered <- df_bubble %>%
#   group_by(Species, DisplayLabel, salinity_class) %>%
#   summarise(MeanRA = mean(RA), .groups = "drop") %>%
#   arrange(salinity_class, desc(MeanRA)) %>%
#   pull(DisplayLabel)
# 
# df_bubble <- df_bubble %>%
#   mutate(DisplayLabel = factor(DisplayLabel, levels = rev(unique(sp_ordered))))
# 
# p2 <- ggplot(df_bubble,
#              aes(x = Zone, y = DisplayLabel,
#                  size = RA, colour = salinity_class)) +
#   geom_point(alpha = 0.82, stroke = 0.25, shape = 16) +
#   facet_wrap(~ YearLabel, nrow = 1, strip.position = "top") +
#   scale_size_area(
#     max_size = 14,
#     name     = "Relative\nAbundance (%)",
#     breaks   = c(1, 5, 10, 25, 50),
#     labels   = function(x) paste0(x, "%"),
#     guide    = guide_legend(
#       override.aes = list(colour = "#444444"),
#       order = 2
#     )
#   ) +
#   scale_colour_manual(
#     values = salinity_colors,
#     name   = "Salinity Tolerance",
#     guide  = guide_legend(
#       override.aes = list(size = 4.5),
#       order = 1
#     )
#   ) +
#   labs(
#     title    = "Taxon-Level Community Composition",
#     subtitle = "Bubble area proportional to relative abundance; colour = Salinity Tolerance",
#     x        = "Zone",
#     y        = NULL
#   ) +
#   theme_npg_pub() +
#   theme(
#     axis.text.y      = element_text(size = BASE_FONT - 2.5, face = "italic",
#                                     lineheight = 0.85),
#     axis.text.x      = element_text(angle = 30, hjust = 1),
#     panel.grid.major = element_line(colour = "#EBEBEB", linewidth = 0.3),
#     panel.grid.minor = element_blank()
#   )
# 
# n_sp  <- length(levels(df_bubble$DisplayLabel))
# p2_h  <- max(8, n_sp * 0.32 + 3)
# p2_w  <- max(12, length(yr_order) * 3.5 + 4)
# 
# out2_pdf <- file.path(OUTPUT_DIR, "Fig2_bubble_plot.pdf")
# out2_png <- file.path(OUTPUT_DIR, "Fig2_bubble_plot.png")
# ggsave(out2_pdf, p2, width = p2_w, height = p2_h, device = cairo_pdf)
# ggsave(out2_png, p2, width = p2_w, height = p2_h, dpi = DPI)
# cat("  ✓ Saved:", out2_pdf, "&", out2_png, "\n")

# # ============================================================
# #  PLOT 3 — ALLUVIAL DIAGRAMS (one per Zone)
# # ============================================================
# cat("Building Plot 3: Alluvial Diagrams ...\n")
# 
# p3_list <- list()
# 
# for (z in levels(df_long$Zone)) {
#   
#   # FIX 1: sum within each sample before averaging, same as Plot 1.
#   # FIX 2: apply MIN_ABUND filter AFTER summing so that small species
#   #         are still counted inside "Other" before any row is dropped.
#   #         In the original code the filter ran first, silently removing
#   #         those reads and pulling totals below 100%.
#   df_al <- df_long %>%
#     filter(Zone == z) %>%
#     group_by(ColName, YearLabel, SpeciesGrp, DisplayLabel) %>%
#     summarise(RA = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%   # ← sum per sample
#     group_by(YearLabel, SpeciesGrp, DisplayLabel) %>%
#     summarise(RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%        # ← then average
# #    filter(RA >= MIN_ABUND) %>%                                         # ← filter after aggregation
# #    filter(RA > 0) %>%
#     complete(YearLabel, SpeciesGrp, fill = list(RA = 0)) %>%
#     # left_join(
#     #   df_long %>% distinct(SpeciesGrp, DisplayLabel),
#     #   by = "SpeciesGrp"
#     # ) %>%
#     # mutate(
#     #   DisplayLabel = factor(
#     #     ifelse(SpeciesGrp == "Other", "Other", label_lookup[SpeciesGrp]),
#     #     levels = c(label_lookup[top_species], "Other")
#     #   )
#     # )
#     mutate(
#       DisplayLabel = factor(
#         ifelse(SpeciesGrp == "Other", "Other", label_lookup[SpeciesGrp]),
#         levels = c(label_lookup[top_species], "Other")
#       )
#     )#%>%
#     # mutate(
#     #   DisplayLabel = coalesce(DisplayLabel.x, DisplayLabel.y),
#     #   DisplayLabel = factor(DisplayLabel, levels = c(label_lookup[top_species], "Other"))
#     # ) %>%
#     #select(-DisplayLabel.x, -DisplayLabel.y)
#   
#   total_per_yr <- df_al %>%
#     group_by(YearLabel) %>%
#     summarise(TotRA = sum(RA), .groups = "drop")
#   
#   df_al <- df_al %>%
#     left_join(total_per_yr, by = "YearLabel") %>%
#     mutate(
#       Pct        = RA / TotRA * 100,
#       StratLabel = if_else(Pct >= 3, as.character(DisplayLabel), "")
#     )
#   
#   p3 <- ggplot(df_al,
#                aes(x        = YearLabel,
#                    y        = RA,
#                    alluvium = SpeciesGrp,
#                    stratum  = SpeciesGrp,
#                    fill     = SpeciesGrp)) +
#     geom_flow(
#       stat          = "alluvium",
#       lode.guidance = "frontback",
#       alpha         = 0.55,
#       colour        = NA,
#       aes(fill = SpeciesGrp)
#     ) +
#     geom_stratum(
#       width     = 0.7,
#       colour    = "white",
#       linewidth = 0.3
#     ) +
# #    geom_label(
# #      stat   = "stratum",
# #      aes(label = after_stat(stratum)),
# #      size   = 2.2,
# #      fill   = alpha("white", 0.65),
# #      colour = "#1A1A1A",
# #      label.size    = 0,
# #      label.padding = unit(1.2, "pt"),
# #      fontface      = "italic",
# #      check_overlap = TRUE
# #    ) +
#     scale_fill_species() +
#     scale_y_continuous(
#       labels = function(x) paste0(round(x, 0), "%"),
#       expand = expansion(mult = c(0, 0.02))
#     ) +
#     labs(
#       title    = paste0("Community Dynamics: ", z),
#       # subtitle = paste0(
#       #   "Alluvial flow of relative abundance across sampling Years"
#       #   #, if (MIN_ABUND > 0) paste0("  |  Taxa < ", MIN_ABUND, "% hidden") else ""
#       # ),
#       x = "",
#       y = "Relative Abundance (%)"
#     ) +
#     theme_npg_pub() +
#     theme(
#       axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.85),
#       panel.grid  = element_blank(),
#       #legend.position = "none"
#     )
#   
#   p3_list[[z]] <- p3
#   
#   out3_base <- file.path(OUTPUT_DIR, paste0("Fig3_alluvial_zone_", make.names(z)))
#   ggsave(paste0(out3_base, ".pdf"), p3,
#          width = max(9, length(yr_order) * 2.2 + 2),
#          height = 7, device = cairo_pdf)
#   ggsave(paste0(out3_base, ".png"), p3,
#          width = max(9, length(yr_order) * 2.2 + 2),
#          height = 7, dpi = DPI)
#   cat("  ✓ Zone", z, "saved.\n")
# }



# ============================================================
#  PLOT 3 — ALLUVIAL DIAGRAMS (one per Zone)
# ============================================================
cat("Building Plot 3: Alluvial Diagrams ...\n")

p3_list <- list()

# 1. Define the exact order of zones and their corresponding alphabet prefixes
target_zones <- c("Zone_1", "Zone_2A", "Zone_3", "Zone_2B")
prefixes     <- c("a) ", "b) ", "c) ", "d) ")

# 2. Loop by index instead of factor levels
for (i in seq_along(target_zones)) {
  
  z      <- target_zones[i]
  pfx    <- prefixes[i]
  p_title <- paste0(pfx, z) # Creates "a) Zone_1", "b) Zone_2A", etc.
  
  df_al <- df_long %>%
    filter(Zone == z) %>%
    group_by(ColName, YearLabel, SpeciesGrp, DisplayLabel) %>%
    summarise(RA = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%   
    group_by(YearLabel, SpeciesGrp, DisplayLabel) %>%
    summarise(RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%        
    complete(YearLabel, SpeciesGrp, fill = list(RA = 0)) %>%
    mutate(
      DisplayLabel = factor(
        ifelse(SpeciesGrp == "Other", "Other", label_lookup[SpeciesGrp]),
        levels = c(label_lookup[top_species], "Other")
      )
    ) 
  
  total_per_yr <- df_al %>%
    group_by(YearLabel) %>%
    summarise(TotRA = sum(RA), .groups = "drop")
  
  df_al <- df_al %>%
    left_join(total_per_yr, by = "YearLabel") %>%
    mutate(
      Pct        = RA / TotRA * 100,
      StratLabel = if_else(Pct >= 3, as.character(DisplayLabel), "")
    )
  
  p3 <- ggplot(df_al,
               aes(x        = YearLabel,
                   y        = RA,
                   alluvium = SpeciesGrp,
                   stratum  = SpeciesGrp,
                   fill     = SpeciesGrp)) +
    geom_flow(
      stat          = "alluvium",
      lode.guidance = "frontback",
      alpha         = 0.55,
      colour        = NA,
      aes(fill = SpeciesGrp)
    ) +
    geom_stratum(
      width     = 0.70, # Kept your preferred bar width
      colour    = "white",
      linewidth = 0.3
    ) +
    scale_fill_species() +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 0), "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title    = p_title, # ← UPDATED DYNAMIC TITLE HERE
      # subtitle = "Alluvial flow of relative abundance across sampling Years",
      x        = "",
      y        = ""
    ) +
    theme_npg_pub() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.85),
      panel.grid  = element_blank()
    )
  
  p3_list[[z]] <- p3
  
  # out3_base <- file.path(OUTPUT_DIR, paste0("Fig3_alluvial_zone_", make.names(z)))
  # ggsave(paste0(out3_base, ".pdf"), p3, width = max(9, length(yr_order) * 2.2 + 2), height = 7, device = cairo_pdf)
  # ggsave(paste0(out3_base, ".png"), p3, width = max(9, length(yr_order) * 2.2 + 2), height = 7, dpi = DPI)
  # cat("  ✓ Zone", z, "saved.\n")
}

# ── Optional: combined alluvial panel ────────────────────────
#if (length(p3_list) > 1) {
#  p3_combined <- wrap_plots(p3_list, ncol = 1, guides = "collect") &
#    theme(legend.position = "right")
#  
#  out3c_pdf <- file.path(OUTPUT_DIR, "Fig3_alluvial_ALL_zones.pdf")
#  out3c_png <- file.path(OUTPUT_DIR, "Fig3_alluvial_ALL_zones.png")
#  ggsave(out3c_pdf, p3_combined,
#         width  = max(9, length(yr_order) * 2.2 + 2),
#         height = 7 * length(p3_list),
#         device = cairo_pdf)
#  ggsave(out3c_png, p3_combined,
#         width  = max(9, length(yr_order) * 2.2 + 2),
#         height = 7 * length(p3_list),
#         dpi    = DPI, limitsize = FALSE)
#  cat("  ✓ Combined alluvial panel saved.\n")
#}


# if (length(p3_list) > 1) {
#   p3_combined <- wrap_plots(p3_list, ncol = 1, guides = "collect") &
#     theme(legend.position = "right")   # ← one legend for all panels
#   
#   out3c_pdf <- file.path(OUTPUT_DIR, "Fig3_alluvial_ALL_zones.pdf")
#   out3c_png <- file.path(OUTPUT_DIR, "Fig3_alluvial_ALL_zones.png")
#   ggsave(out3c_pdf, p3_combined,
#          width  = max(9, length(yr_order) * 2.2 + 2),
#          height = 7 * length(p3_list),
#          device = cairo_pdf)
#   ggsave(out3c_png, p3_combined,
#          width  = max(9, length(yr_order) * 2.2 + 2),
#          height = 7 * length(p3_list),
#          dpi    = DPI, limitsize = FALSE)
#   cat("  ✓ Combined alluvial panel saved.\n")
# }

# ── Combined alluvial panel (2x2 Grid Layout) ────────────────────────
if (length(p3_list) > 1) {
  
  p3_combined <- wrap_plots(p3_list, ncol = 2, guides = "collect") &
    theme(legend.position = "right")   # Keep one shared legend on the right
  
  out3c_pdf <- file.path(OUTPUT_DIR, "Fig3_alluvial_2x2_zones.pdf")
  out3c_png <- file.path(OUTPUT_DIR, "Fig3_alluvial_2x2_zones.png")
  
  # 2. Adjust saving dimensions for the grid
  # Individual panel width is calculated based on years; we double it for 2 columns.
  
  # Slightly narrowed the width matrix to force columns closer together
  single_width <- max(5, length(yr_order) * 1 + 2)
  grid_width   = (single_width * 2)
  grid_height  = 6.5 * 2
  
  
  # single_width <- max(9, length(yr_order) * 2.2 + 2)
  # grid_width   <- (single_width * 2) - 1  # Subtracted 1 to account for shared legend space efficiency
  # grid_height  <- 7 * 2                    # 2 rows high instead of 4 rows high
  
  ggsave(out3c_pdf, p3_combined,
         width  = grid_width,
         height = grid_height,
         device = cairo_pdf)
  
  ggsave(out3c_png, p3_combined,
         width  = grid_width,
         height = grid_height,
         dpi    = DPI, limitsize = FALSE)
  
  cat("  ✓ Combined 2x2 alluvial panel saved.\n")
}

# ── DONE ─────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════\n")
cat(" All figures saved to:", normalizePath(OUTPUT_DIR), "\n")
cat(" Files:\n")
cat("  Fig1_stacked_bar.{pdf,png}\n")
cat("  Fig2_bubble_plot.{pdf,png}\n")
cat("  Fig3_alluvial_zone_<Z>.{pdf,png}   (one per zone)\n")
cat("  Fig3_alluvial_ALL_zones.{pdf,png}  (combined panel)\n")
cat("══════════════════════════════════════════════════════════\n")