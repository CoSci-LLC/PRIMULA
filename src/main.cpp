// ============================================================================
// Copyright (C) model++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

#include <Config.hpp>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
   Config config(argc, argv);

   Primula model = config.configModel();

   model.ReadSoilDataset("../tests/malonno/Pedologia_25k_MALONNO.csv", "../tests/malonno/RootReinforcement.csv");
   model.GenerateSoilProperties();
   model.CalculateSafetyFactor();

   model.pr_failure_.writeToFile(
      (fs::path(config.outputPath) / fs::path("prob_failure")).string() + config.extension());

   return EXIT_SUCCESS;
}
