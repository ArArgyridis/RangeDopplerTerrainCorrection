# RangeDopplerTerrainCorrection
An algorithm for Sentinel 1A/1B geo-registration

In order to be build gdal, libxml2, and boost should be installed on the system.

Usage:  ./RangeDopplerTerrainCorrection sentinel_metatada.dim dem_file EPSG pixelSize output_directory

dem_file should have the same EPSG as the output image (EPSG parameter)
