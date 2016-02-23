////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Feb 18 20:44:50 EST 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <fstream>


int MakeFakeBy (std::string const OutFileName)
{
  // Open output file
  std::ofstream f(OutFileName.c_str());
  if (!f.is_open()) {
    std::cerr << "ERROR: Cannot open file: " << OutFileName << std::endl;;
  }


  float const ZStart = -0.9;
  float const ZStop  =  0.9;

  int const NPoints = 1000;

  float const StepSize = (ZStop - ZStart) / (NPoints - 1);

  float Z, Bx, By, Bz;
  for (int i = 0; i != NPoints; ++i) {
    Z  = ZStart + StepSize * i;

    Bx = 0;
    By = 0;
    Bz = 0;

    if (Z < 0.0) {
      Bx = 0;
      By = 0;
    } else if (Z < 0.1) {
      Bx = 0;
      By = 1;
    } else if (Z < 0.2) {
      Bx =  1;
      By = -1;
    } else if (Z < 0.3) {
      Bx = -1;
      By = -1;
    } else if (Z < 0.4) {
      Bx = -1;
      By =  1;
    } else if (Z < 0.5) {
      Bx = 1;
      By = 0;
    } else if (Z < 0.6) {
      Bx = 0;
      By = 0;
    } else if (Z < 0.7) {
      Bx =  0;
      By =  0;
    } else if (Z < 0.8) {
      Bx = 0;
      By = 0;
    } else {
    }

    f << Z << " " << Bx << " " << By << " " << Bz << std::endl;
  }




  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " [OutFileName]" << std::endl;
    return 1;
  }

  MakeFakeBy(argv[1]);

  return 0;
}
