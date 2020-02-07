# -*- coding: utf-8 -*-

# Sophie Dufour-Beauséjour
# s.dufour.beausejour@gmail.com
# Doctorante en sciences de l"eau
# INRS-ETE
# Québec

# Computes pixel statistics from an image file over an
# area defined by a shapefile.
# The shapefile may include more than one polygon.
# Images expected in tiff format.
# Requirements: fiona, rasterio, seaborn (if plotting)
# This file is imported in main.py
# Contents:
# - compute_statistics: compute statistics for tiff AOIs in list, write to txt
# - AOIs_from_shp: write tiff for each feature in AOI shapefile

from __future__ import print_function

import os
import fiona
import rasterio
from rasterio.tools.mask import mask
import numpy as np
from natsort import natsorted
import matplotlib.pyplot as plt
import seaborn
from utils import line_pre_adder

def compute_statistics(AOI_paths, results_dir, band_index,
        no_data_value, overwrite=0, plot=1):
    """Compute pixel statistics from specific band (band_index) for
    list of tiff files (AOI_paths), writing to a
    single text file. NaN and no_data_value pixels are
    masked."""

    # Check if stats already computed
    if os.path.exists(os.path.join(results_dir, "stats", "median.txt")):
        print("stats already computed")
        if overwrite:
            print("overwritting...")
        else:
            return

    n = len(AOI_paths)
    stat_names = ["mean","stddev","median","q1","q3","max","min","N"]
    # Initialize arrays for each statistic, in a dictionary
    data = {}
    for stat_name in stat_names:
        data[stat_name] = np.zeros([1,n])
    if plot:
        # prepare fig
        dim = int(1+np.sqrt(n*2))
        fig, axes = plt.subplots(dim,dim,figsize=(20, 20))
        axes = axes.flatten()

    # Walk through AOI_paths
    print("computing stats for "+str(len(AOI_paths))+" AOIs...")
    for i, AOI_path in enumerate(AOI_paths):
        # Read tiff
        dataset = rasterio.open(AOI_path, mode='r+')
        # Read band to an array of pixel values
        image_data = np.ma.array(dataset.read(band_index))
        # Mask the pixels that are NaN
        nan_mask = np.isnan(image_data)
        nan_masked_data = np.ma.array(image_data, mask=nan_mask)
        # Mask pixels that are no_data value
        masked_data = np.ma.array(nan_masked_data, mask=(nan_masked_data == no_data_value))
        if plot:
            # tiff
            res = 1 # res = 10 for lower resolution
            sample_x = np.arange(0,len(masked_data[0,:]),res)
            sample_y = np.arange(0,len(masked_data[:,0]),res)
            sample_X, sample_Y = np.meshgrid(sample_x, sample_y)
            sample_data = masked_data[sample_Y, sample_X]
            axes[2*i].pcolormesh(sample_x/res, sample_y/res,
                    np.flipud(sample_data),cmap="gray",vmin=-22, vmax=-5)
            plt.axes(axes[2*i]).set_aspect('equal', 'datalim')

        # Compress pixels and compute stats
        masked_data = np.ma.compressed(masked_data)
        data["mean"][0,i] = masked_data.mean()
        data["stddev"][0,i] = np.std(masked_data)
        data["q1"][0,i] = np.percentile(masked_data,25)
        data["q3"][0,i] = np.percentile(masked_data,75)
        data["max"][0,i] = np.max(masked_data)
        data["min"][0,i] = np.min(masked_data)
        data["N"][0,i] = len(masked_data)
        try :
            data["median"][0,i] = np.median(masked_data)
        except :
            "User warning: invalid value encountered in median."
        if plot:
            # histogram
            seaborn.distplot(masked_data,ax=axes[2*i+1],color="k")
            axes[2*i].set_xlabel(str(i+1)+" (dB)")

    # Save each statistic to text file combining all AOIs
    for key in data.keys():
        header = np.hstack(["# pixel "+key+" for list of AOI_tiffs",
            "\n#", map(str, range(n+1))[1:]])
        # Write to text file
        save_path = os.path.join(results_dir, "stats", key+".txt")
        np.savetxt(save_path, data[key], fmt='%.8f', delimiter='\t')
        # Add header
        line_pre_adder(save_path," \t   ".join(header))
    if plot:
        # write plot
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        save_path = os.path.join(results_dir, "stats", "figure.png")
        plt.savefig(save_path, fmt="png", dpi=200)

def AOIs_from_shp(image_path, shapefile_path, results_dir, overwrite=0):
    """From a tiff image (path: image_path), write a tiff (in results_dir)
    for each feature in a shapefile (path: shapefile_path)."""

    # Check if files already there
    if os.path.exists(os.path.join(results_dir, "tiff_AOIs", "AOI_1.tif")):
        print("tiff_AOIs already written")
        if overwrite:
            print("overwritting...")
        else:
            return

    # Read shapefile
    with fiona.open(shapefile_path, "r") as shapefile:
        # Walk through shapefile features
        for i,feature in enumerate(shapefile):
            AOI_path = os.path.join(results_dir,"tiff_AOIs", "AOI_"+str(i+1)+".tif")
            # Get feature geometry
            geoms = [feature["geometry"]]
            # Open image, mask with geometry
            with rasterio.open(image_path) as src:
                out_image, out_transform = mask(src, geoms, crop=True)
                out_meta = src.meta.copy()
            # Write the image (masked and cropped to the geometry)
            out_meta.update({"driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform})
            with rasterio.open(AOI_path, "w", **out_meta) as dest:
                dest.write(out_image)
