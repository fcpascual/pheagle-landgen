## 01_mantel.R - individual-based isolation by distance (spatial subset)

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "R"
source(file.path(script_dir, "config.R"))

need("dartRverse")
suppressPackageStartupMessages(library(dartRverse))
od <- out_dir("mantel")
gl_sp <- read_gl("pe_gl_spatial")
msg("IBD on", nInd(gl_sp), "individuals,", nLoc(gl_sp), "loci")

set.seed(SEED)
pdf(file.path(od, "ibd_plots.pdf"), width = 8, height = 6)
res <- gl.ibd(gl_sp, distance = "euclidean", coordinates = COORD_ARG,
              permutations = N_PERM)
dev.off()

writeLines(capture.output(print(res$mantel)), file.path(od, "mantel_summary.txt"))
m <- res$mantel
write.csv(data.frame(n_ind = nInd(gl_sp), mantel_r = m$statistic,
                     p_value = m$signif, permutations = m$permutations),
          file.path(od, "mantel_summary.csv"), row.names = FALSE)
saveRDS(res, file.path(od, "gl_ibd.rds"))

## extra scatter of genetic vs geographic distance (best effort)
tryCatch({
  pdf(file.path(od, "dgen_vs_dgeo.pdf"), width = 6, height = 6)
  x <- as.vector(res$Dgeo); y <- as.vector(res$Dgen)
  plot(x, y, pch = 19, col = "#3182bd80",
       xlab = "Geographic distance", ylab = "Genetic distance")
  abline(lm(y ~ x), lty = 2)
  dev.off()
}, error = function(e) {
  while (dev.cur() > 2) dev.off()
  msg("Skipped distance scatter:", conditionMessage(e))
})
msg("01_mantel done")
