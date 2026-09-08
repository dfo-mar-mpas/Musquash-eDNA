# ============================================================================
# MULTI-MARKER SPATIAL ANALYSIS: ROBUST ATOMIC VERSION
# ============================================================================

library(raster)
library(sf)
library(gdistance)
library(vegan)
library(adespatial)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. CONFIGURATION
markers <- list(
  #"MiFish" = "E:/projects/Musquash_data/Analysis/results/matrices_input_24/Musq-24-mifish-bray-curtis.tsv",
  #"COI"    = "E:/projects/Musquash_data/Analysis/results/matrices_input_24/Musq-24-COI-bray-curtis.tsv"
  "MiFish" = "E:/projects/Musquash_data/Analysis/results2/COST_PATH/matrices/Musq-25-mifish-aitchison.tsv",
  "COI"    = "E:/projects/Musquash_data/Analysis/results2/COST_PATH/matrices/Musq-25-COI-aitchison.tsv"
)

META_PATH <- "E:/projects/Musquash_data/Analysis/metadata/site_coordinates_metadata.csv"
DEM_PATH  <- "E:/projects/Musquash_data/github repos/musquash_mpa-master/data/Bathymetry/Musquash_dem.tif"

N_ITER      <- 250   
MANTEL_PERM <- 999   
set.seed(123)

# 2. HELPER FUNCTIONS
get_water_dist <- function(sites_df, dem_raster) {
  sites_sf <- st_as_sf(sites_df, coords = c("longitude", "latitude"), crs = 4326)
  sites_utm <- st_transform(sites_sf, crs = st_crs(dem_raster))
  coords_utm <- st_coordinates(sites_utm)
  
  site_depths <- raster::extract(dem_raster, coords_utm)
  thresh <- min(site_depths, na.rm = TRUE) - 0.5
  
  cost_surf <- dem_raster
  values(cost_surf)[values(cost_surf) < thresh] <- NA
  values(cost_surf)[values(cost_surf) >= thresh] <- 1
  
  trans <- transition(cost_surf, transitionFunction = mean, directions = 8)
  trans_corr <- geoCorrection(trans, type = "c")
  
  dist_mat <- costDistance(trans_corr, coords_utm) / 1000
  
  if(any(is.infinite(dist_mat))) {
    euc <- as.matrix(dist(coords_utm) / 1000)
    final_m <- as.matrix(dist_mat)
    final_m[is.infinite(final_m)] <- euc[is.infinite(final_m)] * 1.5
    dist_mat <- as.dist(final_m)
  }
  
  attr(dist_mat, "Labels") <- sites_df$site_id
  return(dist_mat)
}

# 3. ANALYSIS LOOP
all_summaries <- list()

for (m_name in names(markers)) {
  message("\n>>> Analyzing Marker: ", m_name)
  
  beta_in <- read.table(markers[[m_name]], row.names = 1, header = TRUE, sep = "\t")
  meta_in <- read.csv(META_PATH)
  
  common_sites <- intersect(rownames(beta_in), meta_in$site_id)
  beta_mat <- as.matrix(beta_in[common_sites, common_sites])
  meta_df <- meta_in %>% filter(site_id %in% common_sites) %>% arrange(match(site_id, common_sites))
  
  geo_dist <- get_water_dist(meta_df, raster(DEM_PATH))
  geo_mat  <- as.matrix(geo_dist)
  
  sample_sizes <- seq(5, length(common_sites), by = 1)
  marker_iters <- list()
  
  for (k in sample_sizes) {
    message("  Calculating Power for N = ", k)
    # Initialize with numeric zeros to force atomic types
    k_results <- data.frame(
      n_sites = rep(k, N_ITER), 
      mantel_r = numeric(N_ITER), 
      p_val = numeric(N_ITER), 
      adj_r2 = numeric(N_ITER)
    )
    
    for (i in 1:N_ITER) {
      sel <- sample(common_sites, k)
      sub_beta <- as.dist(beta_mat[sel, sel])
      sub_geo  <- as.dist(geo_mat[sel, sel])
      
      suppressMessages({
        # Mantel Test - Force result to numeric
        m_test <- tryCatch({
          vegan::mantel(sub_beta, sub_geo, permutations = MANTEL_PERM)
        }, error = function(e) return(list(statistic = NA, signif = NA)))
        
        k_results$mantel_r[i] <- as.numeric(m_test$statistic)
        k_results$p_val[i]    <- as.numeric(m_test$signif)
        
        # dbRDA - Force result to numeric
        n_axes <- min(k - 2, 3)
        space_pcoa <- cmdscale(sub_geo, k = n_axes)
        db_mod <- tryCatch({
          capscale(sub_beta ~ space_pcoa)
        }, error = function(e) return(NULL))
        
        if(!is.null(db_mod)){
          k_results$adj_r2[i] <- as.numeric(RsquareAdj(db_mod)$adj.r.squared)
        } else {
          k_results$adj_r2[i] <- NA
        }
      })
    }
    marker_iters[[as.character(k)]] <- k_results
  }
  all_summaries[[m_name]] <- do.call(rbind, marker_iters) %>% mutate(Marker = m_name)
}

# 4. AGGREGATE RESULTS
final_results <- do.call(rbind, all_summaries)

# Ensure the columns are strictly numeric before summarizing
final_results$mantel_r <- as.numeric(final_results$mantel_r)
final_results$adj_r2   <- as.numeric(final_results$adj_r2)

summary_table <- final_results %>%
  group_by(Marker, n_sites) %>%
  summarise(
    mean_r = mean(mantel_r, na.rm=T), 
    sd_r = sd(mantel_r, na.rm=T),
    power = mean(p_val < 0.05, na.rm=T),
    mean_r2 = mean(adj_r2, na.rm=T), 
    sd_r2 = sd(adj_r2, na.rm=T),
    .groups = "drop"
  )

# --- 5. IDENTIFY OPTIMAL SAMPLING SIZES ---
cat("\n=== SAMPLING SUFFICIENCY THRESHOLDS (Power >= 0.8) ===\n")
thresholds <- summary_table %>%
  filter(power >= 0.8) %>%
  group_by(Marker) %>%
  slice(1) %>%
  select(Marker, n_sites, power)

print(thresholds)

# --- 6. MASTER PLOT ---
theme_q1 <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), axis.title = element_text(face = "italic"),
          legend.position = "top", strip.background = element_blank())
}

# Use scale_color_brewer for colorblind friendly Q1 palettes
p1 <- ggplot(summary_table, aes(n_sites, mean_r, color=Marker, fill=Marker)) +
  geom_ribbon(aes(ymin=mean_r-sd_r, ymax=mean_r+sd_r), alpha=0.15, color=NA) +
  geom_line(linewidth=1) + 
  labs(subtitle="A. Mantel Correlation (r)", x="No. of Sites", y="r") + 
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  theme_q1()

p2 <- ggplot(summary_table, aes(n_sites, power, color=Marker)) +
  geom_line(linewidth=1) + geom_hline(yintercept=0.8, linetype="dashed", color="grey40") +
  labs(subtitle="B. Detection Power", x="No. of Sites", y="Proportion p < 0.05") + 
  scale_color_brewer(palette="Set1") +
  theme_q1()

p3 <- ggplot(summary_table, aes(n_sites, mean_r2, color=Marker, fill=Marker)) +
  geom_ribbon(aes(ymin=mean_r2-sd_r2, ymax=mean_r2+sd_r2), alpha=0.15, color=NA) +
  geom_line(linewidth=1) + 
  labs(subtitle="C. Explained Variance", x="No. of Sites", y=expression(Adj.~italic(R)^2)) + 
  scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
  theme_q1()

master_plot <- (p1 | p2 | p3) + plot_layout(guides = "collect") +
  plot_annotation(title = "Marker Comparison: Spatial Sufficiency for Musquash Estuary")

#print(master_plot)

# SAVE
#ggsave("E:/projects/Musquash_data/Analysis/results2/COST_PATH/plots/Figure_Marker_Sufficiency_Final.tiff", master_plot, width = 12, height = 4.5, dpi = 600, compression = "lzw")




library(ggsci) # Required for the NPG palette



# --- 6. MASTER PLOT (NPG STYLE) ---
theme_q1 <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), 
          axis.title = element_text(face = "italic"),
          legend.position = "top", 
          strip.background = element_blank(),
          plot.title = element_text(face = "bold", size = 14))
}

p1 <- ggplot(summary_table, aes(n_sites, mean_r, color=Marker, fill=Marker)) +
  geom_ribbon(aes(ymin=mean_r-sd_r, ymax=mean_r+sd_r), alpha=0.15, color=NA) +
  geom_line(linewidth=1) + 
  labs(subtitle="A. Mantel Correlation (r)", x="No. of Sites", y="r") + 
  scale_color_npg() + scale_fill_npg() +
  theme_q1()

p2 <- ggplot(summary_table, aes(n_sites, power, color=Marker)) +
  geom_line(linewidth=1) + 
  geom_hline(yintercept=0.8, linetype="dashed", color="grey40") +
  labs(subtitle="B. Detection Power", x="No. of Sites", y="Proportion p < 0.05") + 
  scale_color_npg() +
  guides(color = "none") +  # This removes the legend for this plot
  theme_q1()

p3 <- ggplot(summary_table, aes(n_sites, mean_r2, color=Marker, fill=Marker)) +
  geom_ribbon(aes(ymin=mean_r2-sd_r2, ymax=mean_r2+sd_r2), alpha=0.15, color=NA) +
  geom_line(linewidth=1) + 
  labs(subtitle="C. Explained Variance", x="No. of Sites", y=expression(Adj.~italic(R)^2)) + 
  scale_color_npg() + scale_fill_npg() +
  theme_q1()

master_plot <- (p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom") #+
  #plot_annotation(title = "Marker Comparison: Spatial Sufficiency for Musquash Estuary")

print(master_plot)

# SAVE
ggsave("E:/projects/Musquash_data/Analysis/results2/COST_PATH/plots/Costpath_25_aitchison.tiff", master_plot, width = 12, height = 4.5, dpi = 600, compression = "lzw")
# ggsave("E:/projects/Musquash_data/Analysis/results2/COST_PATH/plots/CostPath_25_jaccard.pdf", master_plot, width = 12, height = 4.5, dpi = 600)
