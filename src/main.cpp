// ============================================================================
// Copyright (C) model++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

#include <Config.hpp>

int main(int argc, char **argv)
{
   Config config(argc, argv);

   Primula model = config.configModel();


   for (size_t r = 0; r < model.probslope_.nRows; r++) {
      for (size_t c = 0; c < model.probslope_.nCols; c++) {
         if (model.slope_(r, c) == model.slope_.nodata_value)
            model.probslope_(r, c) = model.probslope_.nodata_value;
         else if (model.probslope_(r, c) < 0)
            model.probslope_(r, c) = 1.1028656e-06;
      }
   }

   model.ReadCSV("../tests/malonno/Pedologia_25k_MALONNO.csv", 100);
   model.GenerateLandslides("../tests/malonno/RootReinforcement.csv", 100);
   model.pr_failure_.writeToFile("prob_failure.asc");

   return EXIT_SUCCESS;
}
