# Run this script to generate all tutorial figures.
# Figures are saved to the figures/ folder in the website repo.
#
# data_path should point at a folder containing both IBA_CO1_homogenate_2019_SE.zip
# and IBA_CO1_lysate_2019_SE.zip. Lysate is downsampled to homogenate's size
# right after loading, and every figure after that (map, PCoA, Shannon,
# barplots) uses that single balanced, pooled sample set - one load/merge
# pass total, kept small and fast throughout.

library(asvoccur)
library(vegan)
library(ape)
library(rworldmap)
library(lubridate)

data_path  <- "~/Downloads/IBA/lys_hom"
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

# Filter to sufficient depth and complete dates, then downsample lysate down
# to homogenate's size so both extraction methods are represented equally and
# PCoA's eigendecomposition (O(n^3) in sample count) stays fast. ix is the
# fixed, pooled sample set used by every figure from here on.
sample_depth <- Matrix::colSums(merged$counts)
print(summary(sample_depth))
ix_all     <- which(sample_depth >= 100000 & !is.na(yday))
dataset_id <- factor(sub(":.*", "", colnames(merged$counts)))

set.seed(1)
hom_ix <- ix_all[dataset_id[ix_all] == "IBA_CO1_homogenate_2019_SE"]
lys_ix <- ix_all[dataset_id[ix_all] == "IBA_CO1_lysate_2019_SE"]
ix     <- c(hom_ix, sample(lys_ix, size = length(hom_ix)))
message(length(ix), " samples after filtering and downsampling (", length(hom_ix), " each)")

newmap <- rworldmap::getMap(resolution = "low")

# Many samples share the same/nearby site (repeated visits); jitter spreads
# them into a visible cluster instead of hiding them behind one point
set.seed(1)
lon_jit <- lon + rnorm(length(lon), 0, 0.03)
lat_jit <- lat + rnorm(length(lat), 0, 0.03)

png(fig("3_map.png"), width = 700, height = 700)
par(mar = c(2, 2, 3, 1))
plot(newmap, xlim = c(10, 25), ylim = c(55, 70), asp = 1, main = "Sampling sites")
points(lon_jit[ix], lat_jit[ix], pch = 21, col = "black", bg = color_habitat[ix], cex = 1.8)
legend("bottomleft", bty = "o", bg = "white", legend = names(habitat_cols),
       pt.bg = habitat_cols, pch = 21, pt.cex = 1.6, cex = 1.0,
       title = "Habitat")
dev.off()

# ── Tutorial 4: PCoA (Spearman-based dissimilarity) ───────────────────────────

# The IBA data include project-specific cluster assignments in
# associatedSequences, alongside more general sequence groupings such as
# BOLD BINs. Read from loaded, not merged: merge_data() drops that column
# when combining datasets
cluster_lookup <- unique(data.table::rbindlist(lapply(loaded$asvs, function(x)
  x[, .(taxonID, associatedSequences)])))
cluster_id <- cluster_lookup$associatedSequences[match(rownames(merged$counts), cluster_lookup$taxonID)]

# Check that every ASV maps to one cluster
stopifnot(sum(is.na(cluster_id)) == 0, anyDuplicated(cluster_lookup$taxonID) == 0)

clevels <- unique(cluster_id)
cidx    <- match(cluster_id, clevels)
G <- Matrix::sparseMatrix(i = cidx, j = seq_along(cluster_id), x = 1,
                           dims = c(length(clevels), length(cluster_id)))
cluster_counts <- G %*% merged$counts
rownames(cluster_counts) <- clevels
message(nrow(cluster_counts), " clusters (from ", nrow(merged$counts), " ASVs)")

# Subset to ix, then remove clusters absent from all of these samples - such
# all-zero rows contain no information for comparing the selected samples
counts_filt <- cluster_counts[, ix, drop = FALSE]
counts_filt <- counts_filt[Matrix::rowSums(counts_filt) > 0, , drop = FALSE]
counts_filt <- as.matrix(counts_filt)

# Spearman correlation measures similarity in rank order of cluster
# abundances between samples; convert to a 0-1 dissimilarity
spear_cor  <- cor(counts_filt, method = "spearman")
spear_dist <- (1 - spear_cor) / 2

# Cailliez correction since the dissimilarities are not necessarily Euclidean
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

# Latitude shown at 2-degree intervals; sampling date by month
lat_ticks <- seq(ceiling(min(lat_range) / 2) * 2, max(lat_range), by = 2)

month_ticks  <- sort(tapply(yday[ix], month[ix], min))
month_labels <- month.abb[as.integer(names(month_ticks))]

method      <- dataset_id[ix]
method_cols <- c(IBA_CO1_homogenate_2019_SE = "#e76f51", IBA_CO1_lysate_2019_SE = "#2a9d8f")
method_labs <- c(IBA_CO1_homogenate_2019_SE = "homogenate", IBA_CO1_lysate_2019_SE = "lysate")

xlab <- paste0("PC1 (", round(pcoa_res$values$Rel_corr_eig[1] * 100), "%)")
ylab <- paste0("PC2 (", round(pcoa_res$values$Rel_corr_eig[2] * 100), "%)")

png(fig("4_pcoa.png"), width = 900, height = 900)
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

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     col = "black", bg = method_cols[as.character(method)], pch = 21, cex = 1,
     xlab = xlab, ylab = ylab, main = "Colour by extraction method")
legend("bottomleft", bty = "n", pch = 19, cex = 1, inset = c(1, 0),
       col = method_cols, legend = method_labs)
dev.off()

# ── Tutorial 5: Alpha diversity (Shannon index by habitat) ───────────────────
# Reuses counts_filt from above instead of building a second large dense matrix

shannon <- vegan::diversity(counts_filt, MARGIN = 2)

png(fig("5_shannon.png"), width = 700, height = 500)
par(mar = c(7, 4, 2, 1))
boxplot(shannon ~ habitat[ix], las = 2, xlab = "",
        ylab = "Shannon diversity",
        col = habitat_cols[levels(factor(habitat[ix]))])
dev.off()

# ── Tutorial 6: Barplots ──────────────────────────────────────────────────────
# Identifies the most abundant taxa across the selected samples, calculates
# their mean relative abundance for each habitat and month, and displays the
# results as stacked barplots. Uses the same ix as steps 3-5, so barplots
# reflect the balanced, pooled sample set rather than every sample loaded

plot_barplots <- function(rank, size_taxa = -1, top_x = 10) {
  mycols <- colorRampPalette(c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
    "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
    "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
  ))

  habitats <- names(habitat_cols)
  par(mfrow = c(length(habitats), 1), mar = c(2, 3, 2, 14), xpd = TRUE)

  # Select the most abundant taxa across the samples retained in step 3
  ok <- sort(
    Matrix::rowMeans(cladecounts$norm[[rank]][, ix]),
    index.return = TRUE,
    decreasing = TRUE
  )$ix[1:top_x]

  if (size_taxa == -1) {
    size_taxa <- min(1.5, 8 / length(ok))
  }

  for (hab in habitats) {
    ix_hab <- intersect(ix, which(habitat == hab))

    monthly_averages_matr <- matrix(
      ncol = 12,
      nrow = nrow(cladecounts$norm[[rank]])
    )

    for (j in 1:12) {
      ix2 <- intersect(ix_hab, which(month == j))

      if (length(ix2) > 1) {
        monthly_averages_matr[, j] <-
          Matrix::rowMeans(cladecounts$norm[[rank]][, ix2])
      }

      if (length(ix2) == 1) {
        monthly_averages_matr[, j] <-
          cladecounts$norm[[rank]][, ix2]
      }

      if (length(ix2) == 0) {
        monthly_averages_matr[, j] <- NA
      }
    }

    barplot(
      monthly_averages_matr[ok, ],
      col = mycols(length(ok)),
      main = hab
    )

    legend(
      "bottomleft",
      bty = "n",
      pch = 19,
      col = mycols(length(ok))[length(ok):1],
      cex = size_taxa,
      inset = c(1, 0),
      legend = rownames(cladecounts$norm[[rank]])[rev(ok)]
    )
  }
}

png(fig("6_barplots.png"), width = 800, height = 1400)
plot_barplots(rank = 4, top_x = 8)
dev.off()

message("Done — figures saved to ", figure_dir)
