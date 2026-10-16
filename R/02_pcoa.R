## 02_pcoa.R - PCoA on the full (post-exclusion) set; no coordinates needed

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "R"
source(file.path(script_dir, "config.R"))

need("dartRverse")
suppressPackageStartupMessages(library(dartRverse))
od <- out_dir("pcoa")
pe_gl <- read_gl("pe_gl")

pe_gl <- gl.filter.allna(pe_gl)
msg("after allna filter:", nLoc(pe_gl), "loci")

## Check subsample run
set.seed(1)
pe_gl <- gl.subsample.loci(pe_gl, n = 2000000, verbose = 0)

msg("PCoA on", nInd(pe_gl), "individuals,", nLoc(pe_gl), "loci")

pc <- gl.pcoa(pe_gl, nf = N_PCS, plot.out = FALSE)

eig <- pc$eig
eig_pct <- round(100 * eig / sum(eig), 1)

scores <- as.data.frame(pc$scores[, 1:3])
names(scores) <- c("PC1", "PC2", "PC3")
scores$ind <- indNames(pe_gl)

write.csv(scores, file.path(od, "pcoa_scores.csv"), row.names = FALSE)
write.csv(data.frame(axis = seq_along(eig), eigenvalue = eig, pct = eig_pct),
          file.path(od, "pcoa_eigenvalues.csv"), row.names = FALSE)
saveRDS(pc, file.path(od, "pcoa.rds"))

col1 <- "#3182bd"
pdf(file.path(od, "pcoa_plots.pdf"), width = 12, height = 4.5)

## page 1: three panels
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
pairs_pc <- list(c(1, 2), c(1, 3), c(2, 3))
for (p in pairs_pc) {
  plot(scores[[p[1]]], scores[[p[2]]], col = col1, pch = 19,
       xlab = paste0("PC", p[1], " (", eig_pct[p[1]], "%)"),
       ylab = paste0("PC", p[2], " (", eig_pct[p[2]], "%)"))
  abline(h = 0, v = 0, lty = 3, col = "grey70")
}
par(op)

## page 2: PC1 vs PC2 with individual labels
op <- par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(scores$PC1, scores$PC2, col = col1, pch = 19,
     xlab = paste0("PC1 (", eig_pct[1], "%)"),
     ylab = paste0("PC2 (", eig_pct[2], "%)"),
     main = "PCoA with individual ids")
text(scores$PC1, scores$PC2, labels = scores$ind, pos = 3, cex = 0.6)
abline(h = 0, v = 0, lty = 3, col = "grey70")

## page 3: variance explained
n_show <- min(10, length(eig_pct))
barplot(eig_pct[1:n_show], names.arg = 1:n_show, col = col1,
        xlab = "PCoA axis", ylab = "% variance")
par(op)
dev.off()
msg("02_pcoa done")
