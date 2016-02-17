////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Mon Dec 14 17:36:45 EST 2015
//
////////////////////////////////////////////////////////////////////

#include "TParticleTrajectory.h"

TParticleTrajectory::TParticleTrajectory ()
{
  // Default constructor
}



TParticleTrajectory::~TParticleTrajectory ()
{
  // Destruction!
}



void TParticleTrajectory::Add (double const x, double const y, double const z, double const vx, double const vy, double const vz)
{
  // Add an element to the stored vector
  fPositionVelocity.push_back( std::array<double, 6> { {x, y, z, vx, vy, vz} } );
  return;
}
