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

processing = "vsTSX"
band_indices = list()
overwrite = 0

if processing == "ratios":
    processing_done = "_sub_cal_spk_rat2_TC2"
    band_indices.append([1,"HH"]) # Rasterio starts counting at 1, not 0
    band_indices.append([2,"HV"]) # Rasterio starts counting at 1, not 0
    band_indices.append([3,"VH"]) # Rasterio starts counting at 1, not 0
    band_indices.append([4,"VV"]) # Rasterio starts counting at 1, not 0
    band_indices.append([5,"VVHHRatio"]) # Rasterio starts counting at 1, not 0
    band_indices.append([6,"HVHHRatio"]) # Rasterio starts counting at 1, not 0
    band_indices.append([7,"VHVVRatio"]) # Rasterio starts counting at 1, not 0
elif processing == "HAa":
    processing_done = "_sub_cal_C3_spk_HAa_TC2"
    band_indices.append([1,"Entropy"]) # Rasterio starts counting at 1, not 0
    band_indices.append([2,"Anisotropy"]) # Rasterio starts counting at 1, not 0
    band_indices.append([3,"Alpha"]) # Rasterio starts counting at 1, not 0
elif processing == "vsTSX":
    processing_done = "_sub_cal_spk_rat2_TC2"
    band_indices.append([4,"VV"]) # Rasterio starts counting at 1, not 0

no_data_value = 0
exp_file_name_pattern = re.compile("RS2.*().*")

for bay in ["Salluit","DeceptionBay","Kangiqsujuaq"]:
    if "vsTSX" in processing:
        if not "Deception" in bay:
            continue
        else:
            results_dir = "/Volumes/Carbe/Doctorat/0. Analyse/" \
                          "data_analyzed/image_statistics/RS2paper/"+bay+"/vsTSX"
            shapefile_path = "/Volumes/Carbe/Doctorat/BaieDeception/Donnees/Imagerie/VectorData/special_issue_TSX/"+"AOIs_ordered.shp"
            print("Computed stats for vsTSX with AOIs_ordered.shp")
    else:
        results_dir = "/Volumes/Carbe/Doctorat/0. Analyse/" \
                  "data_analyzed/image_statistics/RS2paper/"+bay
        shapefile_path = "/Volumes/Carbe/Doctorat/BaieDeception/Donnees/Imagerie/VectorData/RS2paper/"+bay[0]+"_AOIs_RS2_WGS84_reproj.shp"
    image_dir = "/Volumes/Carbe/Doctorat/BaieDeception/Donnees/Imagerie/RS2_subAOIs/"+bay
    image_folders = utils.find_matching_file_list(image_dir, exp_file_name_pattern, print_list=0)
    for i, image_folder in enumerate(image_folders):
        print("Image "+str(i+1)+ " out of "+str(len(image_folders))+"...")
        image_path = os.path.join(image_dir, image_folder, image_folder+processing_done+".tif")

        # Get list of paths to tiff_AOIs
        tiff_folder = "tiff_"+processing
        # Write tiff AOIs
        tiffstats.AOIs_from_shp(image_path, shapefile_path, os.path.join(results_dir,tiff_folder))

        AOI_files = natsorted(os.listdir(os.path.join(results_dir, tiff_folder)))
        AOI_files = [x for x in AOI_files if (".tif" in x) and (image_folder in x)]
        AOI_paths = [os.path.join(results_dir,tiff_folder,x) for x in AOI_files]

        for band_index in band_indices:
            if band_index[1] == "HH":
                convert_to_dB = 1
            elif band_index[1] == "HV":
                convert_to_dB = 1
            elif band_index[1] == "VH":
                convert_to_dB = 1
            elif band_index[1] == "VV":
                convert_to_dB = 1
            else:
                convert_to_dB = 0
            tiffstats.compute_statistics(AOI_paths, results_dir, band_index, no_data_value, convert_to_dB=convert_to_dB, overwrite=overwrite)


# Plot shapefile over image
# tiffstats.plot_shp_over_tiff(image_path, shapefile_path, results_dir, band_index, convert_to_dB=1)
# tiffstats.plot_shp_over_tiff(image_path, os.path.join(results_dir, os.path.basename(results_dir)[0]+"_RS2_square_AOIs.shp"), results_dir, band_index, convert_to_dB=1)
