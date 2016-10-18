#include "TSurfacePoints_BBoxSimple.h"





TSurfacePoints_BBoxSimple::TSurfacePoints_BBoxSimple (double const X0, double const Y0, double const Z0, double const LX, double const LY, double const LZ, size_t const NX, size_t const NY, size_t const NZ)
{
  fXStepSize = LX / (double) NX;
  fYStepSize = LY / (double) NY;
  fZStepSize = LZ / (double) NZ;

  fXStart = X0 - LX / 2. + fXStepSize / 2.;
  fYStart = Y0 - LY / 2. + fYStepSize / 2.;
  fZStart = Z0 - LZ / 2. + fZStepSize / 2.;

  fNX = NX;
  fNY = NY;
  fNZ = NZ;

  fLX = LX;
  fLY = LY;
  fLZ = LZ;


  fNPlaneZ = fNX * fNY;
  fNPlaneY = fNX * fNZ;
  fNPlaneX = fNY * fNZ;

  fNPoints = 2 * (fNPlaneX + fNPlaneY + fNPlaneZ);



  fSurfaces[0].Init("XY", fNX, fNY, fLX, fLY, fX0, fY0, fZ0 + fLZ / 2, +1);
  fSurfaces[1].Init("XY", fNX, fNY, fLX, fLY, fX0, fY0, fZ0 - fLZ / 2, -1);
  fSurfaces[2].Init("XZ", fNX, fNZ, fLX, fLZ, fX0, fY0 + fLY / 2, fZ0, +1);
  fSurfaces[3].Init("XZ", fNX, fNZ, fLX, fLZ, fX0, fY0 - fLY / 2, fZ0, -1);
  fSurfaces[4].Init("YZ", fNY, fNZ, fLY, fLZ, fX0 + fLX / 2, fY0, fZ0, +1);
  fSurfaces[5].Init("YZ", fNY, fNZ, fLY, fLZ, fX0 - fLX / 2, fY0, fZ0, -1);

}





TSurfacePoints_BBoxSimple::~TSurfacePoints_BBoxSimple ()
{
  // Destruction!
}






TSurfacePoint const TSurfacePoints_BBoxSimple::GetPoint (size_t const i) const
{

  switch (this->GetSurfaceNumber(i)) {
    case 0:
      return fSurfaces[0].GetPoint(i);
    case 1:
      return fSurfaces[1].GetPoint(i - fNPlaneZ);
    case 2:
      return fSurfaces[2].GetPoint(i - 2 * fNPlaneZ);
    case 3:
      return fSurfaces[3].GetPoint(i - 2 * fNPlaneZ - fNPlaneY);
    case 4:
      return fSurfaces[4].GetPoint(i - 2 * fNPlaneZ - 2 * fNPlaneY);
    case 5:
      return fSurfaces[5].GetPoint(i - 2 * fNPlaneZ - 2 * fNPlaneY - fNPlaneX);
    default:
      throw;
  }



  return TSurfacePoint();
}




size_t TSurfacePoints_BBoxSimple::GetNPoints () const
{
  return fNPoints;
}




size_t TSurfacePoints_BBoxSimple::GetSurfaceNumber (size_t const i) const
{
  // The order shall be +XY -XY +XZ -XZ +YZ -YZ

  if (i < fNPlaneZ) {
    return 0;
  } else if (i < 2 * fNPlaneZ) {
    return 1;
  } else if (i < 2 * fNPlaneZ + fNPlaneY) {
    return 2;
  } else if (i < 2 * fNPlaneZ + 2 * fNPlaneY) {
    return 3;
  } else if (i < 2 * fNPlaneZ + 2 * fNPlaneY + fNPlaneX) {
    return 4;
  } else if (i < 2 * fNPlaneZ + 2 * fNPlaneY + 2 * fNPlaneX) {
    return 5;
  } else {
    throw;
  }
}





double TSurfacePoints_BBoxSimple::GetX1 (size_t const i) const
{
  switch (this->GetSurfaceNumber(i)) {
    case 0:
      return fSurfaces[0].GetX1(i);
    case 1:
      return fSurfaces[1].GetX1(i - fNPlaneZ);
    case 2:
      return fSurfaces[2].GetX1(i - 2 * fNPlaneZ);
    case 3:
      return fSurfaces[3].GetX1(i - 2 * fNPlaneZ - fNPlaneY);
    case 4:
      return fSurfaces[4].GetX1(i - 2 * fNPlaneZ - 2 * fNPlaneY);
    case 5:
      return fSurfaces[5].GetX1(i - 2 * fNPlaneZ - 2 * fNPlaneY - fNPlaneX);
    default:
      throw;
  }
}





double TSurfacePoints_BBoxSimple::GetX2 (size_t const i) const
{
  switch (this->GetSurfaceNumber(i)) {
    case 0:
      return fSurfaces[0].GetX2(i);
    case 1:
      return fSurfaces[1].GetX2(i - fNPlaneZ);
    case 2:
      return fSurfaces[2].GetX2(i - 2 * fNPlaneZ);
    case 3:
      return fSurfaces[3].GetX2(i - 2 * fNPlaneZ - fNPlaneY);
    case 4:
      return fSurfaces[4].GetX2(i - 2 * fNPlaneZ - 2 * fNPlaneY);
    case 5:
      return fSurfaces[5].GetX2(i - 2 * fNPlaneZ - 2 * fNPlaneY - fNPlaneX);
    default:
      throw;
  }
}
