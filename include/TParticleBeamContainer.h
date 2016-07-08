#ifndef GUARD_TParticleBeamContainer_h
#define GUARD_TParticleBeamContainer_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Thu Jun  9 14:46:05 EDT 2016
//
// This class is meant to contain several different particle beams.
// It should be the main class used by the module
//
////////////////////////////////////////////////////////////////////

#include <map>

#include "TVector3D.h"
#include "TParticleA.h"
#include "TParticleBeam.h"
#include "TTwiss.h"


class TParticleBeamContainer
{
  public:
    TParticleBeamContainer ();
    ~TParticleBeamContainer ();

    void AddNewParticleBeam (std::string const& Type, std::string const& Name, TVector3D const& X0, TVector3D const& D0, double const E0, double const T0, double const Current, double const Weight = 1);

    TParticleA GetNewParticle ();

    TParticleBeam& GetParticleBeam(std::string const&);

    size_t GetRandomBeamIndexByWeight () const;

    void Clear ();

  private:

    std::vector<double>        fParticleBeamWeightSums;
    std::vector<TParticleBeam> fParticleBeams;

    std::map<std::string, TParticleBeam> fParticleBeamMap;
};


















#endif
