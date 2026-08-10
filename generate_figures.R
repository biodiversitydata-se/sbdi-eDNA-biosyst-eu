# Run this script to generate all tutorial figures.
# Figures are saved to the figures/ folder in the website repo.
#
# data_path should point at a folder with only IBA_CO1_homogenate_2019_SE.zip -
# this produces the main, colleague-validated figures (map.png, pcoa.png,
# barplots.png). data_path_combined additionally includes
# IBA_CO1_lysate_2019_SE.zip, and is used only for the step 4 bonus figure
# (pcoa_combined.png) via its own independent load_data()/merge_data() call -
# it is kept separate so the main pcoa.png stays fast and matches the
# validated single-dataset result.

library(asvoccur)
library(vegan)
library(ape)
library(rworldmap)
library(lubridate)
library(data.table)

data_path          <- "~/Downloads/IBA/homogenate_only"
data_path_combined <- "~/Downloads/IBA/lys_hom"
figure_dir         <- "~/code/github/sbdi-eDNA-biosyst-eu/figures"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
fig <- function(name) file.path(figure_dir, name)

# ── Tutorial 2: Load and prepare data ────────────────────────────────────────

loaded      <- load_data(data_path)
merged      <- merge_data(loaded)
merged_df   <- convert_to_df(merged)
cladecounts <- sum_by_clade(merged$counts, merged$asvs)

# ── Tutorial 3: Map ───────────────────────────────────────────────────────────

lat   <- as.numeric(merged_df$events$decimalLatitude)
lon   <- as.numeric(merged_df$events$decimalLongitude)

# eventDate is a start/end interval string; a few samples have only a bare
# year (no month/day) - parse_date_time() leaves those as NA rather than
# guessing, so they're excluded downstream instead of getting a misleading date
event_start <- lubridate::parse_date_time(
  sub("/.*", "", merged_df$events$eventDate),
  orders = c("ymd_HMSz", "ymd_HMS", "ymd")
)
month <- lubridate::month(event_start)
yday  <- lubridate::yday(event_start)

# Check ENVO values and adjust habitat_map if needed:
print(unique(merged_df$events$env_local_scale))

habitat_map <- c(
  "forested area [ENVO:00000111]"    = "forest",
  "area of cropland [ENVO:01000892]" = "cropland",
  "wetland area [ENVO:00000043]"     = "wetland",
  "grassland area [ENVO]"            = "grassland",
  "urban biome [ENVO:01000249]"      = "urban",
  "alpine biome [ENVO:01001835]"     = "alpine"
)
habitat <- habitat_map[merged_df$events$env_local_scale]

habitat_cols <- c(
  "forest"    = "#2d6a4f",
  "cropland"  = "#d4a017",
  "wetland"   = "#4895ef",
  "grassland" = "#95d5b2",
  "urban"     = "#adb5bd",
  "alpine"    = "#9b5de5"
)
color_habitat <- habitat_cols[habitat]

newmap <- rworldmap::getMap(resolution = "low")

# Many samples share the same/nearby site (repeated visits); jitter spreads
# them into a visible cluster instead of hiding them behind one point
set.seed(1)
lon_jit <- lon + rnorm(length(lon), 0, 0.03)
lat_jit <- lat + rnorm(length(lat), 0, 0.03)

png(fig("map.png"), width = 700, height = 700)
par(mar = c(2, 2, 3, 1))
plot(newmap, xlim = c(10, 25), ylim = c(55, 70), asp = 1, main = "Sampling sites")
points(lon_jit, lat_jit, pch = 21, col = "black", bg = color_habitat, cex = 1.8)
legend("bottomleft", bty = "n", legend = names(habitat_cols),
       pt.bg = habitat_cols, pch = 21, pt.cex = 1.6, cex = 1.0,
       title = "Habitat")
dev.off()

# ── Tutorial 4: PCoA (Spearman, following iba_microbiome.R) ──────────────────

sample_depth <- Matrix::colSums(merged$counts)
print(summary(sample_depth))
ix <- which(sample_depth >= 100000 & !is.na(yday))
message(length(ix), " of ", length(sample_depth), " samples retained (depth >= 100000, complete date)")

# Necessary adaptation: aggregate ASVs to clusters using the same cluster
# assignments as the original analysis. Read from loaded$asvs, not
# merged$asvs: merge_data() only restricts columns (dropping
# associatedSequences) when it actually merges 2+ datasets - with a single
# dataset, Reduce() never calls the merge function and the column survives
# by accident, so relying on merged$asvs here would break silently later
cluster_lookup <- unique(data.table::rbindlist(lapply(loaded$asvs, function(x)
  x[, .(taxonID, associatedSequences)])))
cluster_id <- cluster_lookup$associatedSequences[match(rownames(merged$counts), cluster_lookup$taxonID)]

clevels <- unique(cluster_id)
cidx    <- match(cluster_id, clevels)
G <- Matrix::sparseMatrix(i = cidx, j = seq_along(cluster_id), x = 1,
                           dims = c(length(clevels), length(cluster_id)))
cluster_counts <- G %*% merged$counts
rownames(cluster_counts) <- clevels
message(nrow(cluster_counts), " clusters (from ", nrow(merged$counts), " ASVs)")

counts_filt <- as.matrix(cluster_counts[, ix])

spear_cor  <- cor(counts_filt, method = "spearman")
spear_dist <- (-1 * spear_cor + 1) / 2

pcoa_res <- pcoa(spear_dist, correction = "cailliez")

ramp <- colorRampPalette(rev(c(
  "#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4"
)))

yday_range <- min(yday[ix], na.rm = TRUE):max(yday[ix], na.rm = TRUE)
color_yday <- ramp(length(yday_range))
yday_index <- yday - min(yday_range) + 1

lat_range <- floor(min(lat[ix], na.rm = TRUE)):ceiling(max(lat[ix], na.rm = TRUE))
color_lat <- ramp(length(lat_range))
lat_index <- round(lat) - min(lat_range) + 1

# Latitude legend matches the original analysis's fixed 2-degree step.
# Year-day uses month names instead of a numeric step: point colour still
# varies continuously by day, but month starts are easier to read and
# adapt automatically to whatever period the data covers
lat_ticks <- seq(ceiling(min(lat_range) / 2) * 2, max(lat_range), by = 2)

month_ticks  <- sort(tapply(yday[ix], month[ix], min))
month_labels <- month.abb[as.integer(names(month_ticks))]

xlab <- paste0("PC1 (", round(pcoa_res$values$Rel_corr_eig[1] * 100), "%)")
ylab <- paste0("PC2 (", round(pcoa_res$values$Rel_corr_eig[2] * 100), "%)")

png(fig("pcoa.png"), width = 900, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 6), xpd = TRUE)

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     col = "black", bg = color_yday[yday_index[ix]], pch = 21, cex = 1,
     xlab = xlab, ylab = ylab, main = "Colour by year-day")
legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
       col = color_yday[month_ticks - min(yday_range) + 1],
       legend = month_labels)

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     col = "black", bg = color_habitat[ix], pch = 21, cex = 1,
     xlab = xlab, ylab = ylab, main = "Colour by habitat")
legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
       col = habitat_cols, legend = names(habitat_cols))

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     col = "black", bg = color_lat[lat_index[ix]], pch = 21, cex = 1,
     xlab = xlab, ylab = ylab, main = "Colour by latitude")
legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
       col = color_lat[lat_ticks - min(lat_range) + 1],
       legend = lat_ticks)

barplot(pcoa_res$values$Rel_corr_eig[1:20],
        ylab = "Variance explained", xlab = "Principal coordinate (PC)")
dev.off()

# ── Tutorial 4 bonus: combining datasets (extraction method) ─────────────────
# Independent pipeline over data_path_combined (homogenate + lysate), kept
# separate from the main pipeline above so pcoa.png stays fast and matches
# the validated single-dataset result. Lysate is downsampled to homogenate's
# size, since it alone has 4,000+ samples - too many for PCoA to handle
# quickly, and enough to dominate the plot if left unbalanced.

if (dir.exists(data_path_combined)) {

  loaded_c    <- load_data(data_path_combined)
  merged_c    <- merge_data(loaded_c)
  merged_df_c <- convert_to_df(merged_c)

  lat_c <- as.numeric(merged_df_c$events$decimalLatitude)
  event_start_c <- lubridate::parse_date_time(
    sub("/.*", "", merged_df_c$events$eventDate),
    orders = c("ymd_HMSz", "ymd_HMS", "ymd")
  )
  month_c <- lubridate::month(event_start_c)
  yday_c  <- lubridate::yday(event_start_c)
  habitat_c <- habitat_map[merged_df_c$events$env_local_scale]

  depth_c    <- Matrix::colSums(merged_c$counts)
  ix_c_all   <- which(depth_c >= 100000 & !is.na(yday_c))
  dataset_id <- factor(sub(":.*", "", colnames(merged_c$counts)))

  set.seed(1)
  hom_ix <- ix_c_all[dataset_id[ix_c_all] == "IBA_CO1_homogenate_2019_SE"]
  lys_ix <- ix_c_all[dataset_id[ix_c_all] == "IBA_CO1_lysate_2019_SE"]
  ix_c   <- c(hom_ix, sample(lys_ix, size = length(hom_ix)))
  message(length(ix_c), " samples after downsampling lysate to match homogenate")

  cluster_lookup_c <- unique(data.table::rbindlist(lapply(loaded_c$asvs, function(x)
    x[, .(taxonID, associatedSequences)])))
  cluster_id_c <- cluster_lookup_c$associatedSequences[match(rownames(merged_c$counts), cluster_lookup_c$taxonID)]
  clevels_c <- unique(cluster_id_c)
  cidx_c    <- match(cluster_id_c, clevels_c)
  G_c <- Matrix::sparseMatrix(i = cidx_c, j = seq_along(cluster_id_c), x = 1,
                               dims = c(length(clevels_c), length(cluster_id_c)))
  cluster_counts_c <- G_c %*% merged_c$counts
  rownames(cluster_counts_c) <- clevels_c

  counts_filt_c <- as.matrix(cluster_counts_c[, ix_c])
  spear_cor_c   <- cor(counts_filt_c, method = "spearman")
  spear_dist_c  <- (-1 * spear_cor_c + 1) / 2
  pcoa_res_c    <- pcoa(spear_dist_c, correction = "cailliez")

  yday_range_c <- min(yday_c[ix_c], na.rm = TRUE):max(yday_c[ix_c], na.rm = TRUE)
  color_yday_c <- ramp(length(yday_range_c))
  yday_index_c <- yday_c - min(yday_range_c) + 1

  lat_range_c <- floor(min(lat_c[ix_c], na.rm = TRUE)):ceiling(max(lat_c[ix_c], na.rm = TRUE))
  color_lat_c <- ramp(length(lat_range_c))
  lat_index_c <- round(lat_c) - min(lat_range_c) + 1
  lat_ticks_c <- seq(ceiling(min(lat_range_c) / 2) * 2, max(lat_range_c), by = 2)

  month_ticks_c  <- sort(tapply(yday_c[ix_c], month_c[ix_c], min))
  month_labels_c <- month.abb[as.integer(names(month_ticks_c))]

  method_c      <- dataset_id[ix_c]
  method_cols   <- c(IBA_CO1_homogenate_2019_SE = "#e76f51", IBA_CO1_lysate_2019_SE = "#2a9d8f")
  method_labs   <- c(IBA_CO1_homogenate_2019_SE = "homogenate", IBA_CO1_lysate_2019_SE = "lysate")
  color_habitat_c <- habitat_cols[habitat_c]

  xlab_c <- paste0("PC1 (", round(pcoa_res_c$values$Rel_corr_eig[1] * 100), "%)")
  ylab_c <- paste0("PC2 (", round(pcoa_res_c$values$Rel_corr_eig[2] * 100), "%)")

  png(fig("pcoa_combined.png"), width = 900, height = 900)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 6), xpd = TRUE)

  plot(pcoa_res_c$vectors[, 1], pcoa_res_c$vectors[, 2],
       col = "black", bg = color_yday_c[yday_index_c[ix_c]], pch = 21, cex = 1,
       xlab = xlab_c, ylab = ylab_c, main = "Colour by year-day")
  legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
         col = color_yday_c[month_ticks_c - min(yday_range_c) + 1], legend = month_labels_c)

  plot(pcoa_res_c$vectors[, 1], pcoa_res_c$vectors[, 2],
       col = "black", bg = color_habitat_c[ix_c], pch = 21, cex = 1,
       xlab = xlab_c, ylab = ylab_c, main = "Colour by habitat")
  legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
         col = habitat_cols, legend = names(habitat_cols))

  plot(pcoa_res_c$vectors[, 1], pcoa_res_c$vectors[, 2],
       col = "black", bg = color_lat_c[lat_index_c[ix_c]], pch = 21, cex = 1,
       xlab = xlab_c, ylab = ylab_c, main = "Colour by latitude")
  legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
         col = color_lat_c[lat_ticks_c - min(lat_range_c) + 1], legend = lat_ticks_c)

  plot(pcoa_res_c$vectors[, 1], pcoa_res_c$vectors[, 2],
       col = "black", bg = method_cols[as.character(method_c)], pch = 21, cex = 1,
       xlab = xlab_c, ylab = ylab_c, main = "Colour by extraction method")
  legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
         col = method_cols, legend = method_labs)
  dev.off()
} else {
  message("Skipping combined-dataset bonus figure: data_path_combined not found")
}

# ── Tutorial 5: Barplots ──────────────────────────────────────────────────────

plot_barplots <- function(rank, size_taxa = -1, top_x = 10) {
  mycols <- colorRampPalette(c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
    "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
    "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
  ))
  habitats <- names(habitat_cols)
  par(mfrow = c(length(habitats), 1), mar = c(2, 3, 2, 14), xpd = TRUE)
  ok <- sort(Matrix::rowMeans(cladecounts$norm[[rank]]),
             index.return = TRUE, decreasing = TRUE)$ix[1:top_x]
  if (size_taxa == -1) { size_taxa <- min(1.5, 8 / length(ok)) }
  for (hab in habitats) {
    ix <- which(habitat == hab)
    monthly_averages_matr <- matrix(ncol = 12,
                                    nrow = nrow(cladecounts$norm[[rank]]))
    for (j in 1:12) {
      ix2 <- intersect(ix, which(month == j))
      if (length(ix2) > 1) monthly_averages_matr[, j] <- Matrix::rowMeans(cladecounts$norm[[rank]][, ix2])
      if (length(ix2) == 1) monthly_averages_matr[, j] <- cladecounts$norm[[rank]][, ix2]
      if (length(ix2) == 0) monthly_averages_matr[, j] <- 0
    }
    barplot(monthly_averages_matr[ok, ], col = mycols(length(ok)), main = hab)
    legend("bottomleft", bty = "n", pch = 19,
           col    = mycols(length(ok))[length(ok):1],
           cex    = size_taxa, inset = c(1, 0),
           legend = rownames(cladecounts$norm[[rank]])[rev(ok)])
  }
}

png(fig("barplots.png"), width = 800, height = 1400)
plot_barplots(rank = 4, top_x = 8)
dev.off()

message("Done — figures saved to ", figure_dir)
