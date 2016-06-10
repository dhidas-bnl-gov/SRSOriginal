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





void TParticleBeamContainer::AddNewParticleBeam (std::string const& Type, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const E0, double const T0, double const Current, double const Weight)
{
  if (fParticleBeamMap.count(Name) != 0) {
    throw;
  }

  if (fParticleBeamWeightSums.size() == 0) {
    fParticleBeamWeightSums.push_back(Weight);
  } else {
    fParticleBeamWeightSums.push_back(fParticleBeamWeightSums.back() + Weight);
  }

  fParticleBeams.push_back( TParticleBeam(Type, X0, D0, E0, T0, Current) );
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
  // UPDATE: needs to actually be implemented
  return 0;
}



