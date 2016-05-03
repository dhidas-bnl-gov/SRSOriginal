#include "TSurfacePoints_Rectangle.h"

TSurfacePoints_Rectangle::TSurfacePoints_Rectangle (int const& nx1, int const& nx2, double const& x1, double const& x2, double const& x, double const& y, double const& z)
{
  // Constructor

  fNX1 = nx1;
  fNX2 = nx2;
  fX1  = x1;
  fX2  = x2;

  fX0 = x;
  fY0 = y;
  fZ0 = z;
}


TSurfacePoints_Rectangle::~TSurfacePoints_Rectangle ()
{
  // Destructor
}




TSurfacePoint const TSurfacePoints_Rectangle::GetPoint (size_t const i) const
{
  TSurfacePoint p;
  return p;
}




size_t TSurfacePoints_Rectangle::GetNPoints () const
{
  return 0;
}
