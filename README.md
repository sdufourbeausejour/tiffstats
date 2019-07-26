# tiffstats

[![DOI](https://zenodo.org/badge/199067255.svg)](https://zenodo.org/badge/latestdoi/199067255)

tiffstats is a Python script for computing pixels statistics from a TIFF image
over shapefile features

## Usage

```python
import tiffstats

tiffstats.tiff_AOIs_from_shp(image_path, shapefile_path, results_dir) # writes a tiff for each feature/AOI in shapefile
tiffstats.compute_statistics(AOI_paths, results_dir, band_index, no_data_value) # computes pixel stats for each AOI tiff
```

## Contributing
Pull requests are welcome.

## License
[MIT](https://choosealicense.com/licenses/mit/)
