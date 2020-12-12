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

   primula.slope_.Read("../tests/malonno/slope_MALONNO.asc");
   primula.slope_.Print("slope.asc", "asc", 8);
   primula.twi_.Read("../tests/malonno/twi_MALONNO.asc");
   primula.twi_.Print("twi.asc", "asc", 8);
   primula.soil_type_.Read("../tests/malonno/soils_MALONNO_v2.asc");
   primula.soil_type_.Print("soils.asc", "asc", 8);
   primula.soil_depth_.Read("../tests/malonno/soils_MALONNO.asc");
   primula.soil_depth_.Print("depth.asc", "asc", 8);
   primula.dusaf_.Read("../tests/malonno/dusaf_MALONNO.asc");
   primula.dusaf_.Print("dusaf.asc", "asc", 8);
   primula.probslope_.Read("../tests/malonno/PROBSLOPE_MALONNO.asc");
   primula.probslope_.Print("probslope.asc", "asc", 8);
   for (size_t i = 0; i < primula.probslope_.attribute_.size(); i++) {
      auto &c = primula.probslope_.attribute_.at(i);
      if (primula.slope_.attribute_.at(i) == primula.slope_.nodata_value_)
         c = primula.probslope_.nodata_value_;
      else if (c < 0)
         c = 1.1028656e-06;
   }

   primula.ReadCSV("../tests/malonno/Pedologia_25k_MALONNO.csv", 100);
   primula.GenerateLandslides("../tests/malonno/RootReinforcement.csv", 100);
   primula.pr_failure_.Print("prob_failure.asc", "asc", 4);

   return EXIT_SUCCESS;
}
