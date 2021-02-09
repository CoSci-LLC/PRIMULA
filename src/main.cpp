/**
 *  Copyright (c) 2020-2021 CoSci LLC, USA <software@cosci-llc.com>
 *  
 *  This file is part of PRIMULA.
 *  
 *  PRIMULA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  PRIMULA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with PRIMULA.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <Config.hpp>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
   Config config(argc, argv);

   Primula model = config.configModel();

   model.GenerateSoilProperties();
   model.CalculateSafetyFactor();

   model.pr_failure_.writeToFile(
      (fs::path(config.outputPath) / fs::path("prob_failure")).string() + config.extension());

   return EXIT_SUCCESS;
}
