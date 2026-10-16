## config.R - shared settings and helpers. Sourced by every step script.
## Every setting can be overridden with an environment variable (PE_*),
## see slurm/run_pe_popgen.sbatch.

env_chr <- function(key, default) {
  v <- Sys.getenv(key, unset = "")
  if (nzchar(v)) v else default
}
env_num <- function(key, default) {
  v <- Sys.getenv(key, unset = "")
  if (nzchar(v)) as.numeric(v) else default
}

## ---- inputs / outputs ----
VCF         <- env_chr("PE_VCF",  "data/philippine_eagle_35.vcf")
META        <- env_chr("PE_META", "data/pe_metadata.csv")   # cols: id, latitude, longitude
OUTDIR      <- env_chr("PE_OUTDIR", "results")
EXCLUDE_IDS <- trimws(strsplit(env_chr("PE_EXCLUDE", "03"), ",")[[1]])

## ---- coordinates ----
## PE_COORD_MODE: "latlon" (decimal degrees) or "xy" (projected, e.g. UTM,
## Euclidean distances in whatever units the table uses).
## PE_XCOL / PE_YCOL: metadata columns used as X (lon/Easting) and Y (lat/Northing).
COORD_MODE <- env_chr("PE_COORD_MODE", "latlon")
X_COL      <- env_chr("PE_XCOL", "longitude")
Y_COL      <- env_chr("PE_YCOL", "latitude")
if (!COORD_MODE %in% c("latlon", "xy")) stop("PE_COORD_MODE must be latlon or xy", call. = FALSE)
COORD_ARG  <- COORD_MODE                       # passed to gl.ibd(coordinates = )
XLAB <- if (COORD_MODE == "xy") "X (Easting)"  else "Longitude"
YLAB <- if (COORD_MODE == "xy") "Y (Northing)" else "Latitude"

## ---- general ----
NCPUS <- as.integer(env_num("SLURM_CPUS_PER_TASK", 1))
SEED  <- as.integer(env_num("PE_SEED", 1))

## ---- optional locus filters applied after exclusions (NA = skip) ----
FILTER_CALLRATE <- env_num("PE_CALLRATE", NA)   # e.g. 0.95
FILTER_MAF      <- env_num("PE_MAF", NA)        # e.g. 0.05

## ---- 01_mantel ----
N_PERM <- as.integer(env_num("PE_NPERM", 999))

## ---- 02_pcoa ----
N_PCS <- as.integer(env_num("PE_NPCS", 5))

## ---- 03_tess3 ----
KMAX          <- as.integer(env_num("PE_KMAX", 6))
TESS_REP      <- as.integer(env_num("PE_TESS_REP", 20))
TESS_MAXIT    <- as.integer(env_num("PE_TESS_MAXIT", 200))
TESS_MAX_SNPS <- as.integer(env_num("PE_TESS_MAX_SNPS", 100000))
K_PLOT <- {
  v <- env_chr("PE_K_PLOT", "")
  if (nzchar(v)) as.integer(strsplit(v, ",")[[1]]) else 2:KMAX
}

## ---- output folder layout ----
DIRS <- c(data = "00_data", mantel = "01_mantel", pcoa = "02_pcoa",
          tess3 = "03_tess3", map = "04_map")

out_dir <- function(key) {
  d <- file.path(OUTDIR, DIRS[[key]])
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## ---- helpers ----
msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n", sep = " ")
}

need <- function(pkgs) {
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop("Missing R packages in this environment: ",
         paste(pkgs[!ok], collapse = ", "), call. = FALSE)
  }
}

read_gl <- function(name) {
  f <- file.path(OUTDIR, DIRS[["data"]], paste0(name, ".rds"))
  if (!file.exists(f)) stop("Missing ", f, ". Run 00_prep first.", call. = FALSE)
  readRDS(f)
}

## Non-interactive: absorb stray plots so no Rplots.pdf is written
if (!interactive()) pdf(NULL)
