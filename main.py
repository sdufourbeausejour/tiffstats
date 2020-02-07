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
import re
from natsort import natsorted
# from this project
import tiffstats
import utils


band_index = [1,"HH"] # Rasterio starts counting at 1, not 0
no_data_value = 0
exp_file_name_pattern = re.compile("RS2.*().*")

for bay in ["Salluit", "DeceptionBay", "Kangiqsujuaq"]:
    results_dir = "/Users/sdufourbeausejour/Desktop/Doctorat/0. Analyse/" \
                  "data_analyzed/image_statistics/RS2paper/"+bay
    shapefile_path = "/Users/sdufourbeausejour/Desktop/Doctorat/0. Analyse/" \
                     "data_analyzed/image_statistics/RS2paper/"+bay+"/"+bay[0]+"_RS2_square_AOIs.shp"
    image_dir = "/Volumes/Crabe/Doctorat/BaieDeception/Donnees/Imagerie/RS2_subAOIs/"+bay

    image_folders = utils.find_matching_file_list(image_dir, exp_file_name_pattern, print_list=0)
    for i, image_folder in enumerate(image_folders):
        print("Image "+str(i+1)+ " out of "+str(len(image_folders))+"...")
        image_path = os.path.join(image_dir, image_folder, image_folder+"_sub_cal_spk_rat_TC.tif")

        # Write tiff AOIs
        tiffstats.AOIs_from_shp(image_path, shapefile_path, results_dir)

        # Get list of paths to tiff_AOIs
        AOI_files = natsorted(os.listdir(os.path.join(results_dir, "tiff_AOIs")))
        AOI_files = [x for x in AOI_files if (".tif" in x) and (image_folder in x)]
        AOI_paths = [os.path.join(results_dir,"tiff_AOIs",x) for x in AOI_files]

        # Compute stats for each AOI in list, save to txt
        tiffstats.compute_statistics(AOI_paths, results_dir, band_index, no_data_value, convert_to_dB=1)

# Plot shapefile over image
# tiffstats.plot_shp_over_tiff(image_path, shapefile_path, results_dir, band_index, convert_to_dB=1)
# tiffstats.plot_shp_over_tiff(image_path, os.path.join(results_dir, os.path.basename(results_dir)[0]+"_RS2_square_AOIs.shp"), results_dir, band_index, convert_to_dB=1)

# # Modify AOI shapefiles
# image_path = "/Volumes/Crabe/Doctorat/BaieDeception/Donnees/Imagerie/" \
#              "RS2_subAOIs/Salluit/RS2_20151219/RS2_20151219_sub_cal_spk_rat_TC.tif"
# shapefile_path = "/Volumes/Crabe/Doctorat/BaieDeception/Donnees/Imagerie/" \
#              "VectorData/RS2_paper/S_AOIs_RS2_WGS84.shp"
# results_dir = "/Users/sdufourbeausejour/Desktop/Doctorat/0. Analyse/" \
#               "data_analyzed/image_statistics/RS2paper/Salluit"
#
# tiffstats.write_new_shp(image_path, shapefile_path, results_dir)

