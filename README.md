# RangeDopplerTerrainCorrection
An algorithm for Sentinel 1A/1B geo-registration

In order to be build gdal, libxml2, and boost should be installed on the system.

Usage:  usage: ./RangeDopplerTerrainCorrection sentinel_metatada.dim dem_file EPSG pixelSize output_file.tif  xMin yMin xMax yMax
               xMin yMin xMax yMax: optional

dem_file should have the same EPSG as the output image (EPSG parameter)
the output file can be written only in .tif format
