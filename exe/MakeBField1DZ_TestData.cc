////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Dec 11 14:56:49 EST 2015
//
// This file is used to generate test data
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <fstream>
#include <cmath>


int MakeBField1DZ_TestData (std::string const OutFileName)
{
  // Constants to change
  float const MaxBy           = 1.0;
  int   const NPoints         = 10000000;
  float const PeriodLength    = 0.021;
  int   const NPeriods        = 61;
  float const UndulatorCenter = 0.0;
  float const PI = 3.14159265359;


  // Calculations
  float const UndulatorLength = PeriodLength * (float) NPeriods;
  float const UndulatorStart  = UndulatorCenter - (UndulatorLength / 2.0);
  float const UndulatorStop   = UndulatorCenter + (UndulatorLength / 2.0);
  float const StepSize = UndulatorLength / (float) (NPoints - 1);

  // Open output file
  std::ofstream f(OutFileName.c_str());


  // Loop in Z and print values to file
  float Z, By, ZFluctuation;
  for (int i = 0; i != NPoints; ++i) {
    ZFluctuation = 2.5 * StepSize * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
    Z  = UndulatorStart + StepSize * (float) i + ZFluctuation;
    By = MaxBy * sin( 2.0 * PI / PeriodLength * Z );
    f << Z << " " << By << "\n";
  }

  // Close file
  f.close();

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName]" << std::endl;
    return 1;
  }

  MakeBField1DZ_TestData(argv[1]);

  return 0;
}
