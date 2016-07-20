#include "TTwiss.h"


TTwiss::TTwiss ()
{
  // Default constructor
}




TTwiss::TTwiss (double const Alpha, double const Beta, double const Epsilon, double const Gamma, TVector3D const& X0)
{
  // Default constructor
}




TTwiss::~TTwiss ()
{
  // Destruction
}



void TTwiss::SetParameters ()
{
  return;
}



void TTwiss::SetX0 (TVector3D const& X0)
{
  fX0 = X0;
  return;
}
