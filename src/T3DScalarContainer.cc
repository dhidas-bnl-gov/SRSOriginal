#include "T3DScalarContainer.h"



T3DScalarContainer::T3DScalarContainer ()
{
  // Default constructor
}




T3DScalarContainer::~T3DScalarContainer ()
{
  // Destruction!!!
}




void T3DScalarContainer::AddPoint (TVector3D const& X, double const V)
{
  fValues.push_back( T3DScalar(X, V) );
  return;
}




size_t T3DScalarContainer::GetNPoints () const
{
  return fValues.size();
}




T3DScalar const& T3DScalarContainer::GetPoint (size_t const i) const
{
  if (i >= fValues.size()) {
    throw;
  }

  return fValues[i];
}







