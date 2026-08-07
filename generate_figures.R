# Run this script to generate all tutorial figures.
# Figures are saved to the figures/ folder in the website repo.
#
# Works unchanged whether data_path points at one COI dataset or several -
# point it at a folder with only IBA_CO1_homogenate_2019_SE.zip first to
# validate against known results, then switch to the combined dataset(s).

library(asvoccur)
library(vegan)
library(ape)
library(rworldmap)
library(lubridate)
library(data.table)

data_path  <- "~/Downloads/IBA"
figure_dir <- "~/code/github/sbdi-eDNA-biosyst-eu/figures"
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
month <- lubridate::month(merged_df$events$eventDate)
yday  <- lubridate::yday(merged_df$events$eventDate)

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
ix <- which(sample_depth >= 100000)
message(length(ix), " of ", length(sample_depth), " samples retained (depth >= 100000)")

# Necessary adaptation: aggregate ASVs to clusters using the same cluster
# assignments as the original analysis (associatedSequences isn't carried
# through merge_data(), so read it from loaded$asvs directly)
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

xlab <- paste0("PC1 (", round(pcoa_res$values$Rel_corr_eig[1] * 100), "%)")
ylab <- paste0("PC2 (", round(pcoa_res$values$Rel_corr_eig[2] * 100), "%)")

png(fig("pcoa.png"), width = 900, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 6), xpd = TRUE)

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     col = "black", bg = color_yday[yday_index[ix]], pch = 21, cex = 1,
     xlab = xlab, ylab = ylab, main = "Colour by year-day")
legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
       col = color_yday[round(seq(1, length(color_yday), length.out = 5))],
       legend = yday_range[round(seq(1, length(yday_range), length.out = 5))])

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     col = "black", bg = color_habitat[ix], pch = 21, cex = 1,
     xlab = xlab, ylab = ylab, main = "Colour by habitat")
legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
       col = habitat_cols, legend = names(habitat_cols))

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     col = "black", bg = color_lat[lat_index[ix]], pch = 21, cex = 1,
     xlab = xlab, ylab = ylab, main = "Colour by latitude")
legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
       col = color_lat[round(seq(1, length(color_lat), length.out = 5))],
       legend = lat_range[round(seq(1, length(lat_range), length.out = 5))])

barplot(pcoa_res$values$Rel_corr_eig[1:20],
        ylab = "Variance explained", xlab = "Principal coordinate (PC)")
dev.off()

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
