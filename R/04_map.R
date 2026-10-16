## 04_map.R - sampling map (static PDF + interactive HTML if leaflet is available)

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(this_file)) dirname(normalizePath(this_file)) else "R"
source(file.path(script_dir, "config.R"))

need("dartRverse")
suppressPackageStartupMessages(library(dartRverse))
od <- out_dir("map")
pe_gl <- read_gl("pe_gl")

pts_df <- data.frame(ind = indNames(pe_gl),
                     lon = pe_gl@other$latlon$lon,
                     lat = pe_gl@other$latlon$lat,
                     stringsAsFactors = FALSE)
pts_df <- pts_df[complete.cases(pts_df[, c("lon", "lat")]), ]
write.csv(pts_df, file.path(od, "map_points.csv"), row.names = FALSE)

pdf(file.path(od, "sampling_map.pdf"), width = 7, height = 7)
plot(pts_df$lon, pts_df$lat, pch = 21, bg = "#3182bd", col = "black", cex = 1.4,
     xlab = XLAB, ylab = YLAB,
     main = paste0("Sampling locations (n = ", nrow(pts_df), ")"))
text(pts_df$lon, pts_df$lat, labels = pts_df$ind, pos = 3, cex = 0.6)
dev.off()

## interactive HTML (not self-contained, so no pandoc needed)
if (COORD_MODE != "latlon") {
  msg("Projected (xy) coordinates; skipped interactive leaflet map")
} else if (requireNamespace("leaflet", quietly = TRUE) &&
    requireNamespace("htmlwidgets", quietly = TRUE)) {
  tryCatch({
    w <- leaflet::leaflet(pts_df)
    w <- leaflet::addTiles(w)
    w <- leaflet::addCircleMarkers(w, lng = ~lon, lat = ~lat, radius = 5,
           stroke = TRUE, weight = 1, color = "black",
           fillColor = "#3182bd", fillOpacity = 0.85,
           popup = ~paste0("<b>", ind, "</b>"))
    htmlwidgets::saveWidget(w, file.path(od, "sampling_map.html"),
                            selfcontained = FALSE)
  }, error = function(e) msg("Skipped leaflet map:", conditionMessage(e)))
} else {
  msg("leaflet/htmlwidgets not installed; skipped interactive map")
}
msg("04_map done")
