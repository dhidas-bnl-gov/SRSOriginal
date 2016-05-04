////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Apr 26 16:58:39 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>

#include "TBFieldUniformB.h"
#include "TBField3DZRegularized.h"
#include "TBFieldSquareWave.h"
#include "TBFieldIdeal1D.h"
#include "TVector3D.h"

int IntegrateField (TBField& B, double const ZBegin, double const ZEnd, double const StepSize = 0.000001)
{
  double I1X = 0;
  double I1Y = 0;
  double I2X = 0;
  double I2Y = 0;
  for (double Z = ZBegin; Z <= ZEnd; Z += StepSize) {
    I1X += B.GetBx(0, 0, Z) * StepSize;
    I1Y += B.GetBy(0, 0, Z) * StepSize;
    I2X += I1X;
    I2Y += I1Y;
    std::cout << Z << " " << I2Y << std::endl;
  }

  //printf("I1X: %12.6E\n", I1X);
  //printf("I1Y: %12.6E\n", I1Y);
  //printf("I2X: %12.6E\n", I2X);
  //printf("I2Y: %12.6E\n", I2Y);

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TBFieldIdeal1D TBF(0.220, 6, 0, 1.000);
  IntegrateField(TBF, -1, 1, 0.001);

  return 0;
}
