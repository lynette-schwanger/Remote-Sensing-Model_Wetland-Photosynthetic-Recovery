# Wetland GPP Modeling with HLS Remote Sensing Data
# Exploratory model: EVI, NDMI, NDWI, VARI
# Final model: EVI, NDMI

library(tidyverse)
library(lubridate)
library(mgcv)
library(terra)

dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/rasters", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)

# read in already cleaned data
train_df <- read_csv("data/train_joined_clean.csv")
test_df  <- read_csv("data/test_joined_clean.csv")

# compute spectral indices
# Band key: B02 = Blue, B03 = Green, B04 = Red, nir = NIR, swir1 = SWIR1

compute_indices <- function(df) {
  df %>%
    mutate(
      EVI  = 2.5 * (nir - B04) / (nir + 6 * B04 - 7.5 * B02 + 1),
      NDMI = (nir - swir1) / (nir + swir1),
      NDWI = (B03 - nir)   / (B03 + nir),
      VARI = (B03 - B04)   / (B03 + B04 - B02))
}

train_df <- compute_indices(train_df)
test_df  <- compute_indices(test_df)

write_csv(train_df, "data/train_final.csv")
write_csv(test_df, "data/test_final.csv")


# fit gam model: EVI + NDMI + NDWI + VARI
gam_exploratory <- gam(
  GPP_NT_VUT_REF ~ s(EVI, k = 30) + s(NDMI, k = 30) +
    s(NDWI, k = 30) + s(VARI, k = 30),
  data   = train_df,
  family = gaussian(),
  method = "REML")

summary(gam_exploratory)

saveRDS(gam_exploratory, "models/gam_evi_ndmi_ndwi_vari.rds")


# observed vs predicted plot for test data
rmse_fun <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

p_test     <- predict(gam_exploratory, newdata = test_df, type = "response")
rmse_test  <- rmse_fun(test_df$GPP_NT_VUT_REF, p_test)

obs  <- test_df$GPP_NT_VUT_REF
pred <- p_test
r2_test <- 1 - sum((obs - pred)^2, na.rm = TRUE) /
  sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)

png("results/figures/obs_vs_pred_test.png", width = 7, height = 7,
    units = "in", res = 300)

lims <- range(c(obs, pred), na.rm = TRUE)
par(mar = c(8, 5, 4, 2) + 0.1)

plot(pred, obs,
     xlim = lims, ylim = lims,
     xlab = expression("Predicted GPP (gC m"^-2*" d"^-1*")"),
     ylab = expression("Observed GPP (gC m"^-2*" d"^-1*")"),
     main = "Observed vs. Predicted GPP: Test Data\nEVI, NDMI, NDWI, VARI",
     pch  = 1)

abline(0, 1, lty = 2)

usr <- par("usr")
text(x = usr[1] + 0.12 * diff(usr[1:2]),
     y = usr[4] - 0.10 * diff(usr[3:4]),
     labels = paste0("RMSE = ", sprintf("%.2f", rmse_test)))
text(x = usr[1] + 0.12 * diff(usr[1:2]),
     y = usr[4] - 0.18 * diff(usr[3:4]),
     labels = paste0("R\u00B2 = ",  sprintf("%.2f", r2_test)))

dev.off()


# variable importance graph
vars   <- c("EVI", "NDMI", "NDWI", "VARI")
p_base <- predict(gam_exploratory, newdata = train_df, type = "response")
rmse_base <- rmse_fun(train_df$GPP_NT_VUT_REF, p_base)

set.seed(33)
B <- 50

imp_rep <- expand.grid(variable = vars, rep = seq_len(B)) |>
  dplyr::as_tibble() |>
  dplyr::rowwise() |>
  dplyr::mutate(
    delta_rmse = {
      dfp          <- train_df
      dfp[[variable]] <- sample(dfp[[variable]])
      p_perm       <- predict(gam_exploratory, newdata = dfp, type = "response")
      rmse_fun(train_df$GPP_NT_VUT_REF, p_perm) - rmse_base
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(variable) |>
  dplyr::summarise(m  = mean(delta_rmse),
                   se = sd(delta_rmse) / sqrt(B),
                   .groups = "drop") |>
  dplyr::arrange(m)

png("results/figures/variable_importance.png", width = 7, height = 5,
    units = "in", res = 300)

par(mar = c(5, 6, 4, 2) + 0.1)
centers <- barplot(imp_rep$m,
                   horiz     = TRUE,
                   names.arg = imp_rep$variable,
                   xlim      = c(0, max(imp_rep$m + imp_rep$se) * 1.15),
                   xlab      = "Increase in RMSE",
                   main      = "Variable Importance")
arrows(imp_rep$m - imp_rep$se, centers,
       imp_rep$m + imp_rep$se, centers,
       lwd = 1.5, angle = 90, code = 3, length = 0.05)

dev.off()


# final model based on variable importance, EVI and NDMI only
gam_final <- gam(
  GPP_NT_VUT_REF ~ s(EVI, k = 30) + s(NDMI, k = 30),
  data   = train_df,
  family = gaussian(),
  method = "REML")

summary(gam_final)
saveRDS(gam_final, "models/gam_evi_ndmi.rds")


# observed vs predicted plot for test data on final model
p_test_final     <- predict(gam_final, newdata = test_df, type = "response")
rmse_test_final  <- rmse_fun(test_df$GPP_NT_VUT_REF, p_test_final)

obs  <- test_df$GPP_NT_VUT_REF
pred_final <- p_test_final
r2_test_final <- 1 - sum((obs - pred_final)^2, na.rm = TRUE) /
  sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)

png("results/figures/obs_vs_pred_test_final.png", width = 7, height = 7,
    units = "in", res = 300)

lims <- range(c(obs, pred_final), na.rm = TRUE)
par(mar = c(8, 5, 4, 2) + 0.1)

plot(pred_final, obs,
     xlim = lims, ylim = lims,
     xlab = expression("Predicted GPP (gC m"^-2*" d"^-1*")"),
     ylab = expression("Observed GPP (gC m"^-2*" d"^-1*")"),
     main = "Observed vs. Predicted GPP: Test Data\nEVI, NDMI",
     pch  = 1)

abline(0, 1, lty = 2)

usr <- par("usr")
text(x = usr[1] + 0.12 * diff(usr[1:2]),
     y = usr[4] - 0.10 * diff(usr[3:4]),
     labels = paste0("RMSE = ", sprintf("%.2f", rmse_test_final)))
text(x = usr[1] + 0.12 * diff(usr[1:2]),
     y = usr[4] - 0.18 * diff(usr[3:4]),
     labels = paste0("R\u00B2 = ",  sprintf("%.2f", r2_test_final)))

dev.off()


# build index raster stacks from raw HLS bands 
build_index_stacks <- function(raw_dir, out_dir, site_short) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # stack individual band tifs by doy 
  all_files <- list.files(raw_dir, pattern = "_aid0001_10N\\.tif$", full.names = TRUE)
  
  doy_id <- sub(".*_doy([0-9]{7})_aid.*", "\\1", basename(all_files))
  files_by_doy <- split(all_files, doy_id)
  
  stacks_by_doy <- lapply(names(files_by_doy), function(d) {
    f <- files_by_doy[[d]]
    b <- sub(".*_(B[0-9A-Z]+)_doy.*", "\\1", basename(f))
    r <- rast(f[order(b)])
    names(r) <- b[order(b)]
    r
  })
  names(stacks_by_doy) <- names(files_by_doy)

  
# rename bands to consistent names
# B8A = NIR for Sentinel-2, B11 = SWIR1, B12 = SWIR2
keep      <- c("B01", "B02", "B03", "B04", "B8A", "B11", "B12")
new_names <- c("B01", "B02", "B03", "B04", "nir", "swir1", "swir2")

# compute indices and save
for (d in names(stacks_by_doy)) {
  r <- stacks_by_doy[[d]]
  
  # check all required bands are present
  if (!all(keep %in% names(r))) {
    message("Skipping DOY ", d, " — missing bands: ",
            paste(setdiff(keep, names(r)), collapse = ", "))
    next
  }
  
  r2 <- r[[keep]]
  names(r2) <- new_names
  
  # compute indices
  evi  <- 2.5 * (r2[["nir"]] - r2[["B04"]]) /
    (r2[["nir"]] + 6 * r2[["B04"]] - 7.5 * r2[["B02"]] + 1)
  ndmi <- (r2[["nir"]] - r2[["swir1"]]) / (r2[["nir"]] + r2[["swir1"]])
  ndwi <- (r2[["B03"]] - r2[["nir"]])   / (r2[["B03"]] + r2[["nir"]])
  vari <- (r2[["B03"]] - r2[["B04"]])   / (r2[["B03"]] + r2[["B04"]] - r2[["B02"]])
  
  # clean up VARI (can produce inf/-inf where denominator = 0)
  vari <- ifel(is.finite(vari), vari, NA)
  
  names(evi)  <- "EVI"
  names(ndmi) <- "NDMI"
  names(ndwi) <- "NDWI"
  names(vari) <- "VARI"
  
  idx <- c(evi, ndmi, ndwi, vari)
  
  out_file <- file.path(
    out_dir,
    paste0("HLSS30_", site_short, "_stack_doy", d, "_EVI_NDMI_NDWI_VARI.tif")
  )
  writeRaster(idx, out_file, overwrite = TRUE)
  cat("Saved:", basename(out_file), "\n")
}
}

# run for both sites
build_index_stacks(
  raw_dir    = "data/raw/HLS_tw1_2018_rasters",
  out_dir    = "data/HLS_tw1_2018_indices",
  site_short = "tw1"
)

build_index_stacks(
  raw_dir    = "data/raw/HLS_tw4_2018_rasters",
  out_dir    = "data/HLS_tw4_2018_indices",
  site_short = "tw4"
)


# predict gpp rasters for 2018 seasonal dates
# export as GeoTIFFs for use in storymap 

# define seasonal dates and labels
seasons <- list(
  list(label = "Winter", date = "2018-01-20", doy = "2018020"),
  list(label = "Spring", date = "2018-04-20", doy = "2018110"),
  list(label = "Summer", date = "2018-08-18", doy = "2018230"),
  list(label = "Fall",   date = "2018-11-06", doy = "2018310")
)

sites <- list(
  list(name = "US-Tw1", short = "tw1"),
  list(name = "US-Tw4", short = "tw4")
)

# fixed color scale across all rasters for comparability
gpp_range <- c(0, 20)

for (site in sites) {
  for (season in seasons) {
    
    raster_path <- sprintf(
      "data/HLS_%s_2018_indices/HLSS30_%s_stack_doy%s_EVI_NDMI_NDWI_VARI.tif",
      site$short, site$short, season$doy
    )
    
    if (!file.exists(raster_path)) {
      message("Raster not found, skipping: ", raster_path)
      next
    }
    
    r_indices <- rast(raster_path)[[c("EVI", "NDMI")]]
    
    # predict GPP
    r_gpp <- terra::predict(r_indices, gam_final, type = "response")
    names(r_gpp) <- "GPP_pred"
    
    mean_gpp <- global(r_gpp, fun = mean, na.rm = TRUE)[[1]]
    cat(site$name, season$label, "— Mean predicted GPP:",
        round(mean_gpp, 2), "gC m-² d-¹\n")
    
    # export GeoTIFF
    tif_out <- sprintf("results/rasters/%s_Predicted_GPP_%s_%s.tif",
                       site$name, season$label, season$date)
    writeRaster(r_gpp, filename = tif_out, overwrite = TRUE)
    
    # export PNG, consistent color scale
    png_out <- sprintf("results/figures/%s_Predicted_GPP_%s_%s.png",
                       site$name, season$label, season$date)
    png(png_out, width = 6, height = 6, units = "in", res = 300)
    plot(r_gpp,
         main   = sprintf("%s %s %s\nMean Predicted GPP = %.2f gC m\u207B\u00B2 d\u207B\u00B9",
                          site$name, season$label, season$date, mean_gpp),
         range  = gpp_range,
         col    = hcl.colors(100, "viridis"))
    dev.off()
  }
}

cat("\nDone. GeoTIFFs in results/rasters/, figures in results/figures/\n")
