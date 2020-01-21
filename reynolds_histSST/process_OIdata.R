## Stephanie Pennington | Created 2019-05-20
## Script to open Reynolds - NOAA Optimum Interpolation (OI) Sea Surface Temperature (SST) V2 (downloaded
## from https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html) and oull out specific sites

library(ncdf4)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

## ----- First, we need to open the data -----
cat("Opening data...")
input_ncfile <- "sst.mnmean.nc" # open Reynolds data
sites <- read_csv("cell_site_lat_longs.csv") # open site locations

# Create data frame with gridbox corner lat/lons for each site
sites_nc <- tibble(site.id = sites$cell_id, lowleft_lat = sites$lat - 0.5, upright_lat = lowleft_lat + 1, 
                      lowleft_lon = sites$long_360 - 0.5, upright_lon = lowleft_lon + 1)

# ## ----- Pull out individual sites and save as separate files with new name -----
# cat("Creating individual netCDF sites...")
# output_basename <- 'oisst_' # this is the beginning portion of each new file
# output_list <- list()
# 
# for(i in 1:nrow(sites_nc)) {
#   output_fullname <- paste0(output_basename, sites_nc$site.id[i], ".nc") # adding unique site name to each file
#   
#   cat("Pulling out site ID:", sites_nc$site.id[i], "...")
#   sellonlatbox('sst', input_ncfile, output_fullname, 
#                sites_nc$lowleft_lon[i], sites_nc$upright_lon[i], sites_nc$lowleft_lat[i], sites_nc$upright_lat[i], nc34 = 4)
# }

cat("Creating individual netCDF sites using CDO...")
output_basename <- 'oisst_'
output_list <- list()

for(i in 1:nrow(sites_nc)) {
    output_fullname <- paste0(output_basename, sites_nc$site.id[i], ".nc") # adding unique site name to each file

    cat("Pulling out site ID:", sites_nc$site.id[i], "...")
    system(paste0("cd ~/Desktop/reynolds_oi; /opt/local/bin/cdo sellonlatbox,", sites_nc$lowleft_lon[i],",",sites_nc$upright_lon[i],",",
                  sites_nc$lowleft_lat[i],",",sites_nc$upright_lat[i]," ",
                  input_ncfile," ",output_fullname))
  }

## ----- Open each site file and format, save as CSV -----
cat("Pulling out variables...")
cat("Saving as CSV...")
ncid <- list.files(path = "../reynolds_oi/", pattern = "^oisst_") %>% lapply(nc_open) # find files, open all

for(j in 1:length(ncid)) {
  sst_filename <- paste0("../reynolds_oi/sst_data/", sites_nc$site.id[j], "_sst.csv") # adding unique site name to each file
  
  # This file parses time as "days since 1800-01-01 00:00:00" so we need to format time
  x <- tibble(YYYYMMDD = as_date(ncvar_get(ncid[[j]], 'time'), origin = "1800-01-01 00:00:00"), 
              SST = ncvar_get(ncid[[j]], 'sst'))
  write.csv(x, sst_filename)
}

cat("All done!")
