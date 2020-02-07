# -*- coding: utf-8 -*-

# Sophie Dufour-Beauséjour
# s.dufour.beausejour@gmail.com
# Université INRS
# Québec, Canada

# Computes pixel statistics from an image file over an
# area defined by a shapefile.
# The shapefile may include more than one polygon.
# Images expected in tiff format.
# Requirements: fiona, rasterio, seaborn (if plotting)
# Other files:
# - tiffstats.py
# - utils.py

from __future__ import print_function
__author__ = "Sophie Dufour-Beauséjour"

import os
from natsort import natsorted
# from this project
import tiffstats

image_path = "/Volumes/Crabe/Doctorat/BaieDeception/Donnees/Imagerie/" \
             "RS2_subAOIs/Salluit/RS2_20170412/RS2_20170412_sub_cal_spk_rat_TC.tif"
shapefile_path = "/Volumes/Crabe/Doctorat/BaieDeception/Donnees/Imagerie/" \
                 "VectorData/RS2_paper/S_AOIs_RS2_WGS84.shp"
results_dir = "local/results/"
band_index = 1 # Rasterio starts counting at 1, not 0
no_data_value = 0

# Write tiff AOIs
tiffstats.AOIs_from_shp(image_path, shapefile_path, results_dir)

# Get list of paths to tiff_AOIs
AOI_files = natsorted(os.listdir(os.path.join(results_dir, "tiff_AOIs")))
AOI_files = [x for x in AOI_files if ".tif" in x]
AOI_paths = [os.path.join(results_dir,"tiff_AOIs",x) for x in AOI_files]

# Compute stats for each AOI in list, save to txt
tiffstats.compute_statistics(AOI_paths, results_dir, band_index, no_data_value)
