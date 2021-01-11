#!/bin/bash

for f in *.asc
do
  filename="${f%.*}"
  gdal_translate $f $filename.tif
  rm -f tmp.*
done
