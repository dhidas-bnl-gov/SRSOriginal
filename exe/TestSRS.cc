////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 17:04:55 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "TSRS.h"


int TestSRS ()
{
  SRS S;

  S.AddMagneticField("epu_linear.dat", "ZBxByBz");
  //S.AddMagneticField("BField_1.5m.dat", "ZBxByBz");


  std::ofstream of("out_By.dat");
  of << std::scientific;

  for (double z = -2.0; z <= 2.0; z += 0.0001) {
    of << z << " " << S.GetBy(0, 0, z) << std::endl;
  }

  of.close();

  S.AddParticleBeam("electron", "Beam_01", 0, 0, -2, 0, 0, 1, 9, 0, 0.500, 1);
  //S.AddParticleBeam("positron", "Beam_02", 0, 0, -2, 0, 0, 1, 6, 0, 0.500, 1);
  std::cout << S.GetParticleBeam("Beam_01") << std::endl;
  //std::cout << S.GetParticleBeam("Beam_02") << std::endl;

  std::cout << S.GetNewParticle() << std::endl;

  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 1) {
    std::cerr << "Usage: " << argv[0] << " " << std::endl;
    return 1;
  }

  TestSRS();

  return 0;
}
