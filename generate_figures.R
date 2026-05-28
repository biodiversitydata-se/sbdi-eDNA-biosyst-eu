# Run this script to generate all tutorial figures.
# Figures are saved to the figures/ folder in the website repo.

install.packages('remotes')
remotes::install_local("~/code/github/asvoccur")

library(asvoccur)
library(ranger)
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

loaded       <- load_data(data_path)
merged       <- merge_data(loaded)
merged_df    <- convert_to_df(merged)
cladecounts <- sum_by_clade(merged$counts, merged$asvs)

# ── Tutorial 3: Map ───────────────────────────────────────────────────────────

DS_malaise    <- grep("homogenate",  rownames(merged_df$events))
DS_soillitter <- grep("soil_litter", rownames(merged_df$events))

lat   <- merged_df$events$decimalLatitude
lon   <- merged_df$events$decimalLongitude
month <- lubridate::month(merged_df$events$eventDate)
yday  <- lubridate::yday(merged_df$events$eventDate)

# Check ENVO values and adjust habitat_map if needed:
print(unique(merged_df$events$env_local_scale))

habitat_map <- c(
  "ENVO:00000091" = "forest",
  "ENVO:00000077" = "cropland",
  "ENVO:00000043" = "wetland",
  "ENVO:00000054" = "grassland",
  "ENVO:01000249" = "urban",
  "ENVO:00000213" = "alpine"
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

pch <- rep(NA, nrow(merged_df$events))
pch[DS_malaise]    <- 21
pch[DS_soillitter] <- 22

newmap <- rworldmap::getMap(resolution = "low")

png(fig("map.png"), width = 800, height = 900)
par(mar = c(2, 2, 2, 2))
plot(newmap, xlim = c(10, 25), ylim = c(55, 70), asp = 1,
     main = "IBA sampling sites")
points(lon, lat, col = "black", bg = color_habitat, pch = pch, cex = 1.2)
legend("bottomleft", bty = "n", legend = names(habitat_cols),
       pt.bg = habitat_cols, pch = 21, pt.cex = 1.4, cex = 0.9,
       title = "Habitat")
legend("bottomright", bty = "n",
       legend = c("Malaise trap", "Soil/litter"),
       pch = c(21, 22), pt.cex = 1.4, cex = 0.9, title = "Method")
dev.off()

# ── Tutorial 4: PCoA ─────────────────────────────────────────────────────────

col_totals  <- Matrix::colSums(merged$counts)
norm_sparse <- Matrix::t(Matrix::t(merged$counts) / col_totals)
norm_counts <- as.matrix(norm_sparse)

bray_dist <- as.matrix(vegdist(t(norm_counts), method = "bray"))
pcoa_res  <- pcoa(bray_dist, correction = "cailliez")

xlab <- paste0("PC1 (", round(pcoa_res$values$Rel_corr_eig[1] * 100), "%)")
ylab <- paste0("PC2 (", round(pcoa_res$values$Rel_corr_eig[2] * 100), "%)")

png(fig("pcoa.png"), width = 1000, height = 500)
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
par(mar = c(5, 5, 2, 1), xpd = TRUE, cex.axis = 1)

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     pch = pch, col = "white", bg = "white", xlab = xlab, ylab = ylab,
     main = "Malaise trap")
points(pcoa_res$vectors[DS_malaise, 1], pcoa_res$vectors[DS_malaise, 2],
       pch = pch[DS_malaise], col = "black",
       bg = color_habitat[DS_malaise], cex = 1.4)

plot(pcoa_res$vectors[, 1], pcoa_res$vectors[, 2],
     pch = pch, col = "white", bg = "white", xlab = xlab, ylab = ylab,
     main = "Soil/litter")
points(pcoa_res$vectors[DS_soillitter, 1], pcoa_res$vectors[DS_soillitter, 2],
       pch = pch[DS_soillitter], col = "black",
       bg = color_habitat[DS_soillitter], cex = 1.4)

plot.new()
legend("center", bty = "n", legend = names(habitat_cols),
       pt.bg = habitat_cols, pch = 21, pt.cex = 1.6, cex = 1.1,
       title = "Habitat")
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
  ok <- sort(rowMeans(cladecounts$norm[[rank]]),
             index.return = TRUE, decreasing = TRUE)$ix[1:top_x]
  if (size_taxa == -1) { size_taxa <- min(1.5, 8 / length(ok)) }
  for (hab in habitats) {
    ix <- which(habitat == hab)
    monthly_averages_matr <- matrix(ncol = 12,
                                    nrow = nrow(cladecounts$norm[[rank]]))
    for (j in 1:12) {
      ix2 <- intersect(ix, which(month == j))
      if (length(ix2) > 1) monthly_averages_matr[, j] <- rowMeans(cladecounts$norm[[rank]][, ix2])
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

# ── Tutorial 6: Random Forest ─────────────────────────────────────────────────

y  <- as.factor(habitat)
ok <- which(!is.na(y))
y  <- y[ok]

X_sparse <- Matrix::t(norm_sparse)[ok, , drop = FALSE]
keep     <- which(Matrix::colSums(X_sparse > 0) / nrow(X_sparse) >= 0.1)
X        <- as.matrix(X_sparse[, keep, drop = FALSE])

set.seed(1)
n        <- nrow(X)
train_ix <- sample(seq_len(n), size = round(0.8 * n))
test_ix  <- setdiff(seq_len(n), train_ix)

X_train <- X[train_ix, , drop = FALSE];  y_train <- y[train_ix]
X_test  <- X[test_ix,  , drop = FALSE];  y_test  <- y[test_ix]

rf   <- ranger(x = X_train, y = y_train, num.trees = 500,
               importance = "permutation")
pred <- predict(rf, data = X_test)$predictions
cat("Accuracy:", round(mean(pred == y_test) * 100), "%\n")

conf      <- table(Predicted = pred, Observed = y_test)
n_classes <- nlevels(y)

png(fig("rf_confusion.png"), width = 650, height = 600)
par(mar = c(9, 9, 3, 2))
image(seq_len(n_classes), seq_len(n_classes), t(conf[n_classes:1, ]),
      axes = FALSE, xlab = "Observed", ylab = "Predicted",
      col  = colorRampPalette(c("white", "#2d6a4f"))(20),
      main = paste0("Accuracy: ", round(mean(pred == y_test) * 100), "%"))
axis(1, at = seq_len(n_classes), labels = colnames(conf), las = 2, cex.axis = 0.9)
axis(2, at = seq_len(n_classes), labels = rev(rownames(conf)), las = 1, cex.axis = 0.9)
for (i in seq_len(nrow(conf)))
  for (j in seq_len(ncol(conf)))
    text(j, n_classes + 1 - i, conf[i, j], cex = 1.1)
dev.off()

imp <- sort(rf$variable.importance, decreasing = TRUE)
png(fig("rf_importance.png"), width = 700, height = 450)
par(mar = c(5, 5, 2, 2))
plot(imp, pch = 16, col = "#2d6a4f",
     xlab = "ASV rank", ylab = "Permutation importance")
dev.off()

# ── Tutorial 7: Air vs. ground ────────────────────────────────────────────────

bray_dist2 <- as.matrix(vegdist(t(cladecounts$norm$genus), method = "bray"))
pcoa_res2  <- pcoa(bray_dist2, correction = "cailliez")

xlab2 <- paste0("PC1 (", round(pcoa_res2$values$Rel_corr_eig[1] * 100), "%)")
ylab2 <- paste0("PC2 (", round(pcoa_res2$values$Rel_corr_eig[2] * 100), "%)")

png(fig("airvsground_pcoa.png"), width = 1000, height = 500)
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
par(mar = c(5, 5, 2, 1), xpd = TRUE, cex.axis = 1)

plot(pcoa_res2$vectors[, 1], pcoa_res2$vectors[, 2],
     pch = pch, col = "white", bg = "white",
     xlab = xlab2, ylab = ylab2, main = "Malaise trap")
points(pcoa_res2$vectors[DS_malaise, 1], pcoa_res2$vectors[DS_malaise, 2],
       pch = pch[DS_malaise], col = "black",
       bg = color_habitat[DS_malaise], cex = 1.4)

plot(pcoa_res2$vectors[, 1], pcoa_res2$vectors[, 2],
     pch = pch, col = "white", bg = "white",
     xlab = xlab2, ylab = ylab2, main = "Soil/litter")
points(pcoa_res2$vectors[DS_soillitter, 1], pcoa_res2$vectors[DS_soillitter, 2],
       pch = pch[DS_soillitter], col = "black",
       bg = color_habitat[DS_soillitter], cex = 1.4)

plot.new()
legend("center", bty = "n", legend = names(habitat_cols),
       pt.bg = habitat_cols, pch = 21, pt.cex = 1.6, cex = 1.1,
       title = "Habitat")
dev.off()

orders          <- rownames(cladecounts$norm$order)
mean_malaise    <- rowMeans(cladecounts$norm$order[, DS_malaise,    drop = FALSE])
mean_soillitter <- rowMeans(cladecounts$norm$order[, DS_soillitter, drop = FALSE])
top             <- order(pmax(mean_malaise, mean_soillitter), decreasing = TRUE)[1:10]

png(fig("airvsground_orders.png"), width = 700, height = 500)
par(mfrow = c(1, 1), mar = c(5, 10, 2, 2))
barplot(rbind(mean_malaise[top], mean_soillitter[top]),
        beside = TRUE, horiz = TRUE, names.arg = orders[top],
        col = c("#4895ef", "#d4a017"), las = 1,
        xlab = "Mean relative abundance",
        main = "Top orders by sampling method")
legend("bottomright", bty = "n",
       legend = c("Malaise trap", "Soil/litter"),
       fill = c("#4895ef", "#d4a017"))
dev.off()

message("Done — figures saved to ", figure_dir)
