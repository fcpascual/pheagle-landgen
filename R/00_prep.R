## 00_prep.R - VCF -> genlight, exclusions, metadata, coordinate-complete subset

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "R"
source(file.path(script_dir, "config.R"))

need("dartRverse")
suppressPackageStartupMessages(library(dartRverse))
od <- out_dir("data")

## ---- read VCF ----
if (!file.exists(VCF)) stop("VCF not found: ", VCF, call. = FALSE)
msg("Reading VCF:", VCF)
pe_gl <- gl.read.vcf(vcffile = VCF)
msg("Read", nInd(pe_gl), "individuals and", nLoc(pe_gl), "loci")

## ---- exclusions ----
drop_ids <- intersect(EXCLUDE_IDS, indNames(pe_gl))
not_found <- setdiff(EXCLUDE_IDS, indNames(pe_gl))
if (length(not_found)) msg("WARNING: exclusion ids not in VCF:", paste(not_found, collapse = ","))
if (length(drop_ids)) {
  pe_gl <- gl.drop.ind(pe_gl, ind.list = drop_ids)
  msg("Dropped:", paste(drop_ids, collapse = ","), "-> n =", nInd(pe_gl))
}

## ---- optional locus filters ----
if (!is.na(FILTER_CALLRATE)) {
  pe_gl <- gl.filter.callrate(pe_gl, method = "loc", threshold = FILTER_CALLRATE,
                              plot.out = FALSE)
  msg("After callrate filter:", nLoc(pe_gl), "loci")
}
if (!is.na(FILTER_MAF)) {
  pe_gl <- gl.filter.maf(pe_gl, threshold = FILTER_MAF, plot.out = FALSE)
  msg("After MAF filter:", nLoc(pe_gl), "loci")
}

## single population for now
pop(pe_gl) <- rep("PE", nInd(pe_gl))

## ---- metadata (ids read as character so leading zeros survive) ----
if (!file.exists(META)) stop("Metadata not found: ", META, call. = FALSE)
meta <- read.csv(META, colClasses = c(id = "character"), stringsAsFactors = FALSE)
need_cols <- c("id", X_COL, Y_COL)
if (!all(need_cols %in% names(meta))) {
  stop("Metadata needs columns: ", paste(need_cols, collapse = ", "), call. = FALSE)
}
if (all(grepl("^[0-9]+$", meta$id))) meta$id <- sprintf("%02d", as.integer(meta$id))
if (anyDuplicated(meta$id)) stop("Duplicate ids in metadata", call. = FALSE)

meta <- meta[match(indNames(pe_gl), meta$id), ]
## X = lon/Easting (X_COL), Y = lat/Northing (Y_COL). In "xy" mode the
## latlon slot simply carries the projected values so later steps can reuse it.
lon <- as.numeric(meta[[X_COL]])
lat <- as.numeric(meta[[Y_COL]])
msg("Coordinate mode:", COORD_MODE, "| X =", X_COL, "| Y =", Y_COL)

pe_gl@other$latlon <- data.frame(lat = lat, lon = lon, row.names = indNames(pe_gl))
pe_gl@other$xy <- data.frame(x = lon, y = lat, row.names = indNames(pe_gl))
pe_gl@other$ind.metrics$lat <- lat
pe_gl@other$ind.metrics$lon <- lon

## ---- coordinate-complete subset ----
has_coords <- complete.cases(pe_gl@other$latlon)
if (sum(has_coords) < 3) stop("Fewer than 3 individuals with coordinates", call. = FALSE)
pe_gl_spatial <- gl.keep.ind(pe_gl, ind.list = indNames(pe_gl)[has_coords])
idx <- match(indNames(pe_gl_spatial), indNames(pe_gl))
pe_gl_spatial@other$latlon <- pe_gl@other$latlon[idx, , drop = FALSE]
pe_gl_spatial@other$xy     <- pe_gl@other$xy[idx, , drop = FALSE]
msg("Full set n =", nInd(pe_gl), "; spatial subset n =", nInd(pe_gl_spatial))

## ---- reports ----
write.csv(data.frame(id = indNames(pe_gl), has_coords = has_coords,
                     lat = lat, lon = lon),
          file.path(od, "sample_table.csv"), row.names = FALSE)

write.csv(as.data.frame(table(chromosome = as.character(chromosome(pe_gl))),
                        stringsAsFactors = FALSE),
          file.path(od, "snps_per_chromosome.csv"), row.names = FALSE)

writeLines(capture.output(print(pe_gl)), file.path(od, "pe_gl_summary.txt"))
writeLines(capture.output(print(pe_gl_spatial)), file.path(od, "pe_gl_spatial_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo.txt"))

saveRDS(pe_gl, file.path(od, "pe_gl.rds"))
saveRDS(pe_gl_spatial, file.path(od, "pe_gl_spatial.rds"))
msg("00_prep done")
