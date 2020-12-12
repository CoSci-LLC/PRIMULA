// ============================================================================
// Copyright (C) PRIMULA++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

#include "primula++.hpp"

int main()
{
   Primula primula;

   primula.slope_ = KiLib::Raster("../tests/malonno/slope_MALONNO.asc");
   primula.slope_.writeToFile("slope.asc");

   primula.twi_ = KiLib::Raster("../tests/malonno/twi_MALONNO.asc");
   primula.twi_.writeToFile("twi.asc");

   primula.soil_type_ = KiLib::Raster("../tests/malonno/soils_MALONNO_v2.asc");
   primula.soil_type_.writeToFile("soils.asc");

   primula.soil_depth_ = KiLib::Raster("../tests/malonno/soils_MALONNO.asc");
   primula.soil_depth_.writeToFile("depth.asc");

   primula.dusaf_ = KiLib::Raster("../tests/malonno/dusaf_MALONNO.asc");
   primula.dusaf_.writeToFile("dusaf.asc");

   primula.probslope_ = KiLib::Raster("../tests/malonno/PROBSLOPE_MALONNO.asc");
   primula.probslope_.writeToFile("probslope.asc");

   for (size_t r = 0; r < primula.probslope_.nRows; r++) {
      for (size_t c = 0; c < primula.probslope_.nCols; c++) {
         if (primula.slope_(r, c) == primula.slope_.nodata_value)
            primula.probslope_(r, c) = primula.probslope_.nodata_value;
         else if (primula.probslope_(r, c) < 0)
            primula.probslope_(r, c) = 1.1028656e-06;
      }
   }

   primula.ReadCSV("../tests/malonno/Pedologia_25k_MALONNO.csv", 100);
   primula.GenerateLandslides("../tests/malonno/RootReinforcement.csv", 100);
   primula.pr_failure_.writeToFile("prob_failure.asc");

   return EXIT_SUCCESS;
}
