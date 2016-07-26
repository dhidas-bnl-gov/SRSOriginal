////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Jun 10 09:21:35 EDT 2016
//             From a coffee shop in Brooklyn, NY
//
////////////////////////////////////////////////////////////////////

#include "TParticleBeamContainer.h"




TParticleBeamContainer::TParticleBeamContainer ()
{
  // Constructor
}





TParticleBeamContainer::~TParticleBeamContainer ()
{
  // Destruction
}





void TParticleBeamContainer::AddNewParticleBeam (std::string const& Type, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const E0, double const T0, double const Current, double const Weight, double const Charge, double const Mass)
{
  if (fParticleBeamMap.count(Name) != 0) {
    throw;
  }

  if (fParticleBeamWeightSums.size() == 0) {
    fParticleBeamWeightSums.push_back(Weight);
  } else {
    fParticleBeamWeightSums.push_back(fParticleBeamWeightSums.back() + Weight);
  }

  if (Type == "custom") {
    fParticleBeams.push_back( TParticleBeam(Type, X0, D0, E0, T0, Current, Weight, Mass) );
  } else {
    fParticleBeams.push_back( TParticleBeam(Type, X0, D0, E0, T0, Current) );
  }
  fParticleBeamMap[Name] = fParticleBeams.back();

  return;
}




TParticleA TParticleBeamContainer::GetNewParticle ()
{
  // UPDATE: incomplete
  return fParticleBeams[ this->GetRandomBeamIndexByWeight() ].GetNewParticle();
}




TParticleBeam& TParticleBeamContainer::GetParticleBeam (std::string const& Name)
{
  // Return a reference to the particle beam given its name
  if (fParticleBeamMap.count(Name) == 0) {
    std::cout << "Didn't find" << std::endl;
    throw;
  }

  return fParticleBeamMap[Name];
}



size_t TParticleBeamContainer::GetRandomBeamIndexByWeight () const
{
  // UPDATE: Get a better random generator with seed setting elsewhere

  // Size of array
  size_t const N = fParticleBeamWeightSums.size();

  // If it's zero we don't really know what we are doing here..
  if (N == 0) {
    throw;
  }

  // If we're 1, that's easy
  if (N == 1) {
    return 0;
  }

  // Get a random double [0, SumOfWeights)
  double const Random = ((double) rand() / (RAND_MAX)) * fParticleBeamWeightSums[N - 1];

  // Not the fastest algorithm, but I guess you don't have thousands of different beams...
  // If you do, let's update this search...
  for (size_t i = 0; i != N; ++i) {
    if (Random < fParticleBeamWeightSums[i]) {
      return i;
    }
  }

  // Just in case you don't find it, something is seriously wrong..
  throw;

  return 0;
}




void TParticleBeamContainer::Clear ()
{
  // Clear the particle beam container contents

  fParticleBeamWeightSums.clear();
  fParticleBeams.clear();
  fParticleBeamMap.clear();

  return;
}
