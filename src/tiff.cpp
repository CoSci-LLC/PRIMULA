#include <libtiff/tiffio.hxx>
#include <spdlog/spdlog.h>
#include <variant>

#include "tiff.hpp"

#define GEOTIFFTAG_MODELPIXELSCALE 33550
#define GEOTIFFTAG_MODELTIEPOINT 33922
#define GEOTIFFTAG_NODATAVALUE 42113

KiLib::Raster fromTiff(const std::string path)
{
        TIFF* tiff = TIFFOpen(path.c_str(), "r");

        if (tiff == NULL)
        {
                spdlog::error("Failed to open {}", path);
                exit(EXIT_FAILURE);
        }

        size_t free_flag = 0;
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
                scaling = new double[2]{1, 1};
                free_flag |= 1;
        }

        if (!TIFFGetField(tiff, GEOTIFFTAG_MODELTIEPOINT, &count, &tiepoint))
        {
                spdlog::critical("Failed to find model tiepoint. Assuming 0, 0");
                tiepoint = new double[6]{0};
                free_flag |= 2;
        }

        if (!TIFFGetField(tiff, GEOTIFFTAG_NODATAVALUE, &count, &nodat))
        {
                spdlog::critical("Failed to find nodata value. Assuming -9999");
                nodat = new char[6];
                free_flag |= 4;
                strcpy(nodat, "-9999\n");
        }

        r.width = (r.nCols - 1) * scaling[0];
        r.height = (r.nRows - 1) * scaling[1];
        r.xllcorner = tiepoint[3];
        r.yllcorner = tiepoint[4] - (r.nRows * scaling[1]);
        r.cellsize = scaling[0];
        r.nodata_value = std::stod(nodat);

        r.data.reserve(r.nCols * r.nRows);

        if (TIFFIsTiled(tiff))
        {
                spdlog::error("There is no support for tiled images yet.");
                exit(EXIT_FAILURE);
        }
        else
        {
                // Rows per strip
                size_t rps = 1;

                // Format is currently undefined: https://www.awaresystems.be/imaging/tiff/tifftags/sampleformat.html
                int16 format = 4;

                // The number of bytes a strip occupies
                uint64 ss = TIFFStripSize64(tiff);

                TIFFGetField(tiff, TIFFTAG_ROWSPERSTRIP, &rps);
                TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &format);

                // Strips are used to avoid having to deal with decompression of data.
                tdata_t buf = _TIFFmalloc(ss);
                for (tstrip_t strip = 0; strip < TIFFNumberOfStrips(tiff); strip += rps)
                {
                        TIFFReadEncodedStrip(tiff, strip, buf, (tsize_t) -1);

                        // There may be issues with in-memory alignment with this approach
                        for (size_t i = 0; i < ss; i += ss/(rps * r.nCols))
                        {
                                switch (format)
                                {
                                case 1:
                                        r.data.emplace_back(*(uint16*)(buf+i));
                                        break;
                                case 2:
                                        r.data.emplace_back(*(int16*)(buf+i));
                                        break;
                                case 3:
                                        r.data.emplace_back(*(double*)(buf+i));
                                        break;
                                default:
                                        spdlog::error("Unknown data format.");
                                        exit(EXIT_FAILURE);
                                }
                        }
                }
                _TIFFfree(buf);
        }

        // Remember to free necessary variables
        if (free_flag & 1)
                delete scaling;

        // Remember to free necessary variables
        if (free_flag & 2)
                delete tiepoint;

        // Remember to free necessary variables
        if (free_flag & 4)
                delete nodat;

        TIFFClose(tiff);

        return r;
}