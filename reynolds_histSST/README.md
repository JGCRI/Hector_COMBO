This read me file describes the directory setup for historic SST calculation from Reynolds NOAA Optimal Interpolation SST Analysis, 
V2 (OISSTV2) 1X1. (Reynolds, R.W., N.A. Rayner, T.M. Smith, D.C. Stokes, and W. Wang, 2002: An improved in situ and satellite SST analysis for climate. J. Climate, 15, 1609-1625)

Extracted SST data is located in `/sst_data/`. History values are monthly and begin Dec-1981 and end Apr-2019 for each site. 

Site location details are located in `cell_site_lat_longs`.

*process_OIdata.R* - opens file sst.mnmean.nc and extracts SST data for each site, saves a raw .nc file to the main directory and a formatted Excel file in /sst_data/.

Contact: Stephanie Pennington stephanie.pennington@pnnl.gov