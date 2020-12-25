#include <libtiff/tiffio.hxx>
#include <spdlog/spdlog.h>
#include <variant>

#include "tiff.hpp"

#define GEOTIFFTAG_MODELPIXELSCALE 33550
#define GEOTIFFTAG_MODELTIEPOINT 33922
#define GEOTIFFTAG_NODATAVALUE 42113

bool isGEOTIFF(const std::string path)
{
        // TODO: implement the tests listed in this section and the ones following it: http://docs.opengeospatial.org/is/19-008r4/19-008r4.html#_tiff_tags_test

        return true;
}

KiLib::Raster fromTiff(const std::string path)
{
        TIFF* tiff = TIFFOpen(path.c_str(), "r");

        if (tiff == NULL)
        {
                spdlog::error("Failed to open {}", path);
                exit(EXIT_FAILURE);
        }

        uint count = 0;
        double* scaling;
        double* tiepoint;
        char* nodat;

        size_t num_dir = 1;

        while (TIFFReadDirectory(tiff))
                num_dir++;

        // The file is can to only contain one directory for now as multiple directories would correspond to multiple rasters
        if (num_dir > 1)
        {
                spdlog::error("Multi-raster images are not supported as of yet.");
                exit(EXIT_FAILURE);
        }

        KiLib::Raster r;

        // Retrieve the width and height of the image. Fail if either can't be retrieved.
        if (!TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &r.nCols) || !TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &r.nRows))
        {
                spdlog::error("Failed to read image width or image height.");
                exit(EXIT_FAILURE);
        }

        if(!TIFFGetField(tiff, GEOTIFFTAG_MODELPIXELSCALE, &count, &scaling))
        {
                spdlog::critical("Failed to find pixel scaling. Assuming 1:1");
                for (size_t i = 0; i < count; i++)
                        scaling[i] = 1;
        }

        if (!TIFFGetField(tiff, GEOTIFFTAG_MODELTIEPOINT, &count, &tiepoint))
        {
                spdlog::critical("Failed to find model tiepoint. Assuming 0, 0");
                for (size_t i = 0; i < count; i++)
                        tiepoint[i] = 0;
        }

        if (!TIFFGetField(tiff, GEOTIFFTAG_NODATAVALUE, &count, &nodat))
        {
                spdlog::critical("Failed to find nodata value. Assuming -9999.");
                nodat = new char[6];
                strcpy(nodat, "-9999\n");
        }

        r.width = r.nCols * scaling[0];
        r.height = r.nRows * scaling[1];
        r.xllcorner = tiepoint[3];
        r.yllcorner = tiepoint[4] - r.height;
        r.cellsize = scaling[0];
        r.nodata_value = std::stod(nodat);

        r.print();

        // TODO: read the data in the tiff file into r.data

        TIFFClose(tiff);

        return r;
}