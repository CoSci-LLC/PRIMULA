// ============================================================================
// Copyright (C) PRIMULA++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++  
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <memory>
#include <string>
#include <algorithm>
#include <vector>
#include <sys/stat.h>
#include <chrono>
#include <limits>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>

//#include "raster.hpp"
//#include "raster_functions.hpp"
//#include "landslide.hpp"
#include "primula++.hpp"

int main(int argc, char **argv)
{
   /*if (argc != 4 && argc != 5) 
   { 
      std::cout << "Usage: slideformap++ <dem file name> <twi file name> [<tree file name>] <number of landslides>" << std::endl;
      return EXIT_SUCCESS;
   }
  
   const std::string demfile = argv[1];
   const std::string twifile = argv[2];
   std::string treefile;
   unsigned int num_landslides;
   if (argc == 5)
   {
      treefile = argv[3];
      num_landslides = std::stoul(argv[4]);
   }
   else
      num_landslides = std::stoul(argv[3]);*/

   Primula primula;
   primula.ReadCSV("../tests/malonno/Pedologia_25k_MALONNO.csv",1000);

   return EXIT_SUCCESS;
}