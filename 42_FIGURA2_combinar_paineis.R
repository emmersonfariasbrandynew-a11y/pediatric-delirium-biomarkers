# 42_FIGURA2_combinar_paineis.R
# See README.md and script_output_map.csv for inputs and outputs.
suppressPackageStartupMessages({ library(magick) })
out_dir <- "output"

a <- image_read(file.path(out_dir, "FIGURA2_PAINELA_REVISED.tif"))
b <- image_read(file.path(out_dir, "FIGURA2_PAINELB_forest_REVISED.tif"))

wa <- image_info(a)$width
wb <- image_info(b)$width
target_w <- max(wa, wb)
a2 <- image_resize(a, paste0(target_w, "x"))
b2 <- image_resize(b, paste0(target_w, "x"))

combined <- image_append(c(a2, b2), stack = TRUE)
image_write(combined, path = file.path(out_dir, "figure2_deliirum_REVISED_CONTEMPORANEO.tif"),
            format = "tiff", compression = "LZW", density = "1200x1200")
cat(sprintf("Combined figure: %dx%d\n", image_info(combined)$width, image_info(combined)$height))
cat("Salvo: output/figure2_deliirum_REVISED_CONTEMPORANEO.tif\n")
