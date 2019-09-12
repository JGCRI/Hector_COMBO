# I want to bring in patterns of temperature change
# multiple by global mean temperature
# put in a format accessiable by the matlab COMBO script

# 9/10/19
# C. Hartin

library(ncdf4)
library(dplyr)
library(tibble)

# testing with GMT
file <- nc_open('C:/Users/hart428/Documents/GitHub/COMBO_Hector/2012.07.09.COMBO/COMBO/Scengen_EPA/PATTERN_tas_ANN_CanESM2_rcp85.nc')
lat <-ncvar_get(file, "lat")
lon <-ncvar_get(file, "lon")
pattern <- ncvar_get(file, "pattern")

#add one spot to longitude values
long <- c(0, lon)

#bring in GMT file for scneario
temp <- read.csv('C:/Users/hart428/Documents/GitHub/COMBO_Hector/2012.07.09.COMBO/COMBO/Scengen_EPA/RCP45/temp_rcp45.csv')

# multiple pattern by GMT for the years needed in COMBO
grid_temp_2000 <- pattern * temp[1,2] 
grid_temp_2000 <- rbind(lat, grid_temp_2000)
grid_temp_2000 <- cbind(long, grid_temp_2000)

grid_temp_2025 <- pattern * temp[26,2]
grid_temp_2025 <- rbind(lat, grid_temp_2025)
grid_temp_2025 <- cbind(long, grid_temp_2025)

grid_temp_2050 <- pattern * temp[51,2]
grid_temp_2050 <- rbind(lat, grid_temp_2050)
grid_temp_2050 <- cbind(long, grid_temp_2050)

grid_temp_2075 <- pattern * temp[76,2]
grid_temp_2070 <- rbind(lat, grid_temp_2075)
grid_temp_2070 <- cbind(long, grid_temp_2075)

grid_temp_2100 <- pattern * temp[101,2]
grid_temp_2100 <- rbind(lat, grid_temp_2100)
grid_temp_2100 <- cbind(long, grid_temp_2100)


