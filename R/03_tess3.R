## 03_tess3.R - TESS3 spatial ancestry on the spatial subset (K = 1..KMAX)
## Needs: tess3r (GitHub: bcm-uga/TESS3_encho_sen), fields

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "R"
source(file.path(script_dir, "config.R"))

need(c("dartRverse", "tess3r", "fields", "sf"))
suppressPackageStartupMessages({ library(dartRverse); library(tess3r) })
od <- out_dir("tess3")
gl_sp <- read_gl("pe_gl_spatial")
ids <- indNames(gl_sp)

## ---- genotype matrix and coordinates ----
gl.matrix <- as.matrix(gl_sp)
rownames(gl.matrix) <- ids

## drop loci that became monomorphic after subsetting
mn <- colMeans(gl.matrix, na.rm = TRUE)
gl.matrix <- gl.matrix[, !is.na(mn) & mn > 0 & mn < 2, drop = FALSE]

## random subsample of loci for speed on large panels
set.seed(SEED)
if (ncol(gl.matrix) > TESS_MAX_SNPS) {
  keep <- sort(sample.int(ncol(gl.matrix), TESS_MAX_SNPS))
  gl.matrix <- gl.matrix[, keep, drop = FALSE]
}

## ---- coordinates ----
## Fit on projected coords (Easting/Northing); convert to lon/lat only for the map.
## In this table, the "lon" column holds Easting and "lat" holds Northing.
coord_xy <- as.matrix(gl_sp@other$latlon[, c("lon", "lat")])
colnames(coord_xy) <- c("x", "y")

if (COORD_MODE == "xy") {
  UTM_EPSG <- as.integer(env_num("PE_EPSG", 32651))   # UTM 51N; use 32650 for zone 50N
  XY_SCALE <- env_num("PE_XY_SCALE", 1000)            # 1000 if table is in km, 1 if meters
  pts <- sf::st_as_sf(data.frame(x = coord_xy[, 1] * XY_SCALE,
                                 y = coord_xy[, 2] * XY_SCALE),
                      coords = c("x", "y"), crs = UTM_EPSG)
  coord_ll <- sf::st_coordinates(sf::st_transform(pts, 4326))
} else {
  coord_ll <- coord_xy
}
colnames(coord_ll) <- c("lon", "lat")
msg("lon range:", paste(round(range(coord_ll[, "lon"]), 2), collapse = " to "),
    "| lat range:", paste(round(range(coord_ll[, "lat"]), 2), collapse = " to "))

## diagnostics
cat("NA fraction per individual:\n"); print(round(rowMeans(is.na(gl.matrix)), 3))
cat("any NaN/Inf:", any(is.nan(gl.matrix)), any(is.infinite(gl.matrix)), "\n")
cat("value range:", range(gl.matrix, na.rm = TRUE), "\n")
cat("coords NA:", anyNA(coord_ll), " duplicated coords:", anyDuplicated(coord_ll), "\n")

## keep loci observed in >= 90% of individuals and with a real minor allele
cr  <- colMeans(!is.na(gl.matrix))
mn  <- colMeans(gl.matrix, na.rm = TRUE)
maf <- pmin(mn, 2 - mn) / 2
gl.matrix <- gl.matrix[, cr >= 0.9 & maf >= 0.05, drop = FALSE]

msg("TESS3 on", nrow(gl.matrix), "individuals and", ncol(gl.matrix), "SNPs")
writeLines(c(paste("individuals:", nrow(gl.matrix)),
             paste("snps_used:", ncol(gl.matrix)),
             paste("seed:", SEED)), file.path(od, "tess3_input_info.txt"))

## ---- run ----
tess_out <- tess3(
  X = gl.matrix,
  coord = coord_xy,
  K = 1:KMAX,
  rep = TESS_REP,
  method = "projected.ls",
  ploidy = 2,
  max.iteration = TESS_MAXIT,
  mask = 0.1,
  keep = "best",
  openMP.core.num = NCPUS
)
saveRDS(tess_out, file.path(od, "tess_out.rds"))

## ---- cross-validation ----
rmse_min <- sapply(1:KMAX, function(k) min(tess_out[[k]]$rmse, na.rm = TRUE))
write.csv(data.frame(K = 1:KMAX, min_cv_rmse = rmse_min),
          file.path(od, "cv_rmse.csv"), row.names = FALSE)

pdf(file.path(od, "tess3_cv_and_barplots.pdf"), width = 8, height = 5)
plot(1:KMAX, rmse_min, type = "b", pch = 19,
     xlab = "K", ylab = "Min CV RMSE", main = "tess3r cross-validation")

Ks <- K_PLOT[K_PLOT >= 1 & K_PLOT <= KMAX]
msg("KMAX =", KMAX, "| plotting K =", paste(Ks, collapse = ","))
for (K in Ks) {
  Q <- qmatrix(tess_out, K = K)
  Qm <- matrix(as.numeric(Q), nrow = nrow(Q), ncol = ncol(Q))
  write.csv(data.frame(ind = ids, Qm, check.names = FALSE),
            file.path(od, sprintf("qmatrix_K%d.csv", K)), row.names = FALSE)

  ord <- order(Qm[, 1], decreasing = TRUE)
  bp <- barplot(t(Qm[ord, , drop = FALSE]), col = hcl.colors(K, "Dark 3"),
                border = NA, space = 0, ylab = "Ancestry proportion",
                main = paste0("tess3r ancestry (K = ", K, ")"),
                names.arg = ids[ord], las = 2, cex.names = 0.6)
}
dev.off()

## ---- ancestry maps (kriging can fail on tiny n; do not kill the job) ----
pdf(file.path(od, "tess3_maps.pdf"), width = 8, height = 7)
for (K in Ks) {
  Q <- qmatrix(tess_out, K = K)
  tryCatch({
    interp.pal <- CreatePalette(hcl.colors(K, "Dark 3"))
    plot(Q, coord_ll, method = "map.max", interpol = FieldsKrigModel(10),
         col = interp.pal, main = paste0("tess3r ancestry map (K = ", K, ")"),
         xlab = XLAB, ylab = YLAB,
         resolution = c(250, 250), cex = 0.5)
  }, error = function(e) {
    plot.new(); text(0.5, 0.5, paste0("Map failed for K = ", K, ":\n", conditionMessage(e)))
    msg("Map failed for K =", K, ":", conditionMessage(e))
  })
}
dev.off()
msg("03_tess3 done")
