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
band_index = [2,"HV"] # Rasterio starts counting at 1, not 0
band_index = [3,"VH"] # Rasterio starts counting at 1, not 0
band_index = [4,"VV"] # Rasterio starts counting at 1, not 0
band_index = [5,"HHVVRatio"] # Rasterio starts counting at 1, not 0
band_index = [6,"HHHVRatio"] # Rasterio starts counting at 1, not 0
band_index = [7,"VVVHRatio"] # Rasterio starts counting at 1, not 0

band_index = [1,"Entropy"] # Rasterio starts counting at 1, not 0
band_index = [2,"Anisotropy"] # Rasterio starts counting at 1, not 0
band_index = [3,"Alpha"] # Rasterio starts counting at 1, not 0


no_data_value = 0
exp_file_name_pattern = re.compile("RS2.*().*")

for bay in ["Salluit","DeceptionBay","Kangiqsujuaq"]:
    results_dir = "/Volumes/Crabe/Doctorat/0. Analyse/" \
                  "data_analyzed/image_statistics/RS2paper/"+bay
    shapefile_path = "/Volumes/Crabe/Doctorat/BaieDeception/Donnees/Imagerie/VectorData/RS2paper/"+bay[0]+"_AOIs_RS2_WGS84_reproj.shp"
    image_dir = "/Volumes/Crabe/Doctorat/BaieDeception/Donnees/Imagerie/RS2_subAOIs/"+bay
    image_folders = utils.find_matching_file_list(image_dir, exp_file_name_pattern, print_list=0)
    for i, image_folder in enumerate(image_folders):
        print("Image "+str(i+1)+ " out of "+str(len(image_folders))+"...")
        image_path = os.path.join(image_dir, image_folder, image_folder+"_sub_cal_C3_spk_HAa_TC2.tif")

        # Write tiff AOIs
        tiffstats.AOIs_from_shp(image_path, shapefile_path, results_dir)

        # Get list of paths to tiff_AOIs
        if "HAa" in image_path:
            tiff_folder = "tiff_AOIs_HAa"
        else:
            tiff_folder = "tiff_AOIs"

        AOI_files = natsorted(os.listdir(os.path.join(results_dir, tiff_folder)))
        AOI_files = [x for x in AOI_files if (".tif" in x) and (image_folder in x)]
        AOI_paths = [os.path.join(results_dir,tiff_folder,x) for x in AOI_files]

        # Compute stats for each AOI in list, save to txt
        if "Ratio" in band_index[1]:
            tiffstats.compute_statistics(AOI_paths, results_dir, band_index, no_data_value, convert_to_dB=0)
        elif "HAa" in image_path:
            tiffstats.compute_statistics(AOI_paths, results_dir, band_index, no_data_value, convert_to_dB=0)
        else:
            tiffstats.compute_statistics(AOI_paths, results_dir, band_index, no_data_value, convert_to_dB=1, overwrite=0)

# Plot shapefile over image
# tiffstats.plot_shp_over_tiff(image_path, shapefile_path, results_dir, band_index, convert_to_dB=1)
# tiffstats.plot_shp_over_tiff(image_path, os.path.join(results_dir, os.path.basename(results_dir)[0]+"_RS2_square_AOIs.shp"), results_dir, band_index, convert_to_dB=1)
