#!/bin/bash

for f in *.tif
do
  filename="${f%.*}"
  gdal_translate -of AAIGrid $f tmp.asc
  gdal_translate -of AAIGrid -co FORCE_CELLSIZE=TRUE tmp.asc $filename.asc
  rm -f tmp.*
done
