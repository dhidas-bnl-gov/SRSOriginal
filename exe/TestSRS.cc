////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 17:04:55 EDT 2016
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>

#include "SRS.h"

#include "TVector3D.h"
#include "TSurfacePoints_RectangleSimple.h"


int TestSRS ()
{
  SRS S;

  S.AddMagneticField("epu_linear.dat", "ZBxByBz");
  //S.AddMagneticField("BField_1.5m.dat", "ZBxByBz");



  S.AddParticleBeam("electron", "Beam_01", TVector3D(0, 0, -1.5), TVector3D(0, 0, 1), 3, 0, 0.500, 1);
  //S.AddParticleBeam("positron", "Beam_02", TVector3D(0, 0, -1.5), TVector3D(0, 0, 1), 3, 0, 0.500, 1);
  std::cout << S.GetParticleBeam("Beam_01") << std::endl;


  TParticleA P = S.GetNewParticle();
  std::cout << P << std::endl;


  // Setup for particle
  S.SetCTStart(0.0);
  S.SetCTStop (3.0);
  S.SetNPointsTrajectory(1901);

  S.SetNewParticle();


  // Calculate the trajectory
  //S.CalculateTrajectory(P);
  //P.GetTrajectory().WriteToFile("del_Trajectory.txt");

  // Calculate the spectrum at a given point
  S.CalculateSpectrum(P, TVector3D(0, 0, 30), 100, 2000, 2000);

  TSurfacePoints_RectangleSimple Surface("XY", 101, 101, 0.080, 0.040, 0, 0, 30, 1);
  //S.CalculatePowerDensity(P, Surface);
  //S.CalculateFlux(P, Surface, 135);




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
