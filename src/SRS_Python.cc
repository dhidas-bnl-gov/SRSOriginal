////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Fri Jun 17 10:39:12 EDT 2016
//
// This is the python-C extension which allows access to the c++
// class SRS.
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "SRS_Python.h"

#include "SRS.h"

#include "TSurfacePoints_RectangleSimple.h"
#include "TSurfacePoints_Rectangle.h"
#include "T3DScalarContainer.h"
#include "TBFieldPythonFunction.h"
#include "TBField3D_Gaussian.h"
#include "TBField3D_Uniform.h"

#include <iostream>
#include <vector>





static void SRS_dealloc(SRSObject* self)
{
  // Python needs to know how to deallocate things in the struct

  delete self->obj;
  self->ob_type->tp_free((PyObject*)self);
}




static PyObject* SRS_new (PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  // Python needs to know how to create things in this struct

  SRSObject* self = (SRSObject*) type->tp_alloc(type, 0);
  if (self != NULL) {

    // Create the new object for self
    self->obj = new SRS();
  }

  // Return myself
  return (PyObject*) self;
}







static TVector3D SRS_ListAsTVector3D (PyObject* List)
{
  TVector3D V;
  if (PyList_Size(List) == 3) {
    Py_INCREF(List);
    V.SetXYZ(PyFloat_AsDouble(PyList_GetItem(List, 0)),
             PyFloat_AsDouble(PyList_GetItem(List, 1)),
             PyFloat_AsDouble(PyList_GetItem(List, 2)));
    Py_DECREF(List);
  } else {
    throw;
  }

  // Return the python list
  return V;
}







static PyObject* SRS_TVector3DAsList (TVector3D const& V)
{
  // Create a python list
  PyObject *PList = PyList_New(0);

  PyList_Append(PList, Py_BuildValue("f", V.GetX()));
  PyList_Append(PList, Py_BuildValue("f", V.GetY()));
  PyList_Append(PList, Py_BuildValue("f", V.GetZ()));

  // Return the python list
  return PList;
}







static PyObject* SRS_Pi (SRSObject* self, PyObject* arg)
{
  return Py_BuildValue("d", TSRS::Pi());
}







static PyObject* SRS_SetCTStart (SRSObject* self, PyObject* arg)
{
  // Set the start time in [m] for caculation

  // Grab the value from input
  double val = PyFloat_AsDouble(arg);

  // Set the object variable
  self->obj->SetCTStart(val);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject* SRS_GetCTStart (SRSObject* self)
{
  // Get the start time in [m] for calculation

  return Py_BuildValue("d", self->obj->GetCTStart());
}




static PyObject* SRS_SetCTStop (SRSObject* self, PyObject* arg)
{
  // Set the stop time in [m] for calculations

  // Grab the value
  double val = PyFloat_AsDouble(arg);

  // Set the object variable
  self->obj->SetCTStop(val);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject* SRS_SetCTStop2 (SRSObject* self, PyObject* args, PyObject *keywds)
{
  // Set the stop time in [m] for calculations

    int voltage;
    char *state = "a stiff";
    char *action = "voom";
    char *type = "Norwegian Blue";

    static char *kwlist[] = {"voltage", "state", "action", "type", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "i|sss", kwlist,
                                     &voltage, &state, &action, &type))
        return NULL;


  self->obj->SetCTStop(0.1);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject* SRS_GetCTStop (SRSObject* self)
{
  // Get the CTStop variable from SRS

  return Py_BuildValue("d", self->obj->GetCTStop());
}




static PyObject* SRS_SetCTStartStop (SRSObject* self, PyObject* args)
{
  // Set the start and stop times for SRS in [m]

  // Grab the values
  double Start, Stop;
  if (! PyArg_ParseTuple(args, "dd", &Start, &Stop)) {
    return NULL;
  }

  // Set the object variable
  self->obj->SetCTStartStop(Start, Stop);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






static PyObject* SRS_GetNPointsTrajectory (SRSObject* self)
{
  // Get the numper of points for trajectory calculaton
  return PyInt_FromSize_t(self->obj->GetNPointsTrajectory());
}




static PyObject* SRS_SetNPointsTrajectory (SRSObject* self, PyObject* arg)
{
  // Set the number of points for trajectory calculation

  // Grab the value from input
  size_t N = PyInt_AsSsize_t(arg);

  // Set the object variable
  self->obj->SetNPointsTrajectory(N);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRS_AddMagneticField (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Set the start and stop times for SRS in [m]
  // UPDATE: Needs comments


  // Grab the values
  char* FileName = "";
  char* FileFormat = "";
  PyObject* LRotation = PyList_New(0);
  PyObject* LTranslation = PyList_New(0);
  PyObject* LScaling = PyList_New(0);

  TVector3D Rotation(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  std::vector<double> Scaling;



  static char *kwlist[] = {"file", "format", "rotations", "translation", "scale", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ss|OOO", kwlist,
                                                           &FileName,
                                                           &FileFormat,
                                                           &LRotation,
                                                           &LTranslation,
                                                           &LScaling)) {
    return NULL;
  }



  if (FileName == "" || FileFormat == "") {
    throw;
  }


  // Check Rotations
  if (PyList_Size(LRotation) != 0) {
    Rotation = SRS_ListAsTVector3D(LRotation);
  }


  // Check Translation
  if (PyList_Size(LTranslation) != 0) {
    Translation = SRS_ListAsTVector3D(LTranslation);
  }


  // Add any scaling factors
  // UPDATE: Check against fileformat number of strings
  for (int i = 0; i < PyList_Size(LScaling); ++i) {
    Scaling.push_back(PyFloat_AsDouble(PyList_GetItem(LScaling, i)));
  }


  // Set the object variable
  self->obj->AddMagneticField(FileName, FileFormat, Rotation, Translation, Scaling);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRS_AddMagneticFieldFunction (SRSObject* self, PyObject* args)
{
  // Set the start and stop times for SRS in [m]

  // Grab the values
  PyObject* Function;
  if (! PyArg_ParseTuple(args, "O:set_callback", &Function)) {
    return NULL;
  }

  Py_INCREF(Function);


  self->obj->AddMagneticField( (TBField*) new TBFieldPythonFunction(Function));

  Py_DECREF(Function);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRS_AddMagneticFieldGaussian (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Set the start and stop times for SRS in [m]


  // Grab the values
  PyObject* LBField = PyList_New(0);
  PyObject* LTranslation = PyList_New(0);
  PyObject* LRotation = PyList_New(0);
  PyObject* LSigma = PyList_New(0);

  TVector3D BField(0, 0, 0);
  TVector3D Sigma(0, 0, 0);
  TVector3D Rotation(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  static char *kwlist[] = {"bfield", "sigma", "rotations", "translation", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|OO", kwlist,
                                                          &LBField,
                                                          &LSigma,
                                                          &LRotation,
                                                          &LTranslation)) {
    return NULL;
  }



  // Check BField
  BField = SRS_ListAsTVector3D(LBField);

  // Check Width
  Sigma = SRS_ListAsTVector3D(LSigma);

  // Check Rotations
  if (PyList_Size(LRotation) != 0) {
    Rotation = SRS_ListAsTVector3D(LRotation);
  }

  // Check Translation
  if (PyList_Size(LTranslation) != 0) {
    Translation = SRS_ListAsTVector3D(LTranslation);
  }


  // Rotate field and sigma
  // UPDATE: check this
  BField.RotateSelfXYZ(Rotation);
  Sigma.RotateSelfXYZ(Rotation);

  // Add field
  self->obj->AddMagneticField( (TBField*) new TBField3D_Gaussian(BField, Translation, Sigma));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRS_AddMagneticFieldUniform (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Set the start and stop times for SRS in [m]
  // UPDATE: Needs comments


  // Grab the values
  PyObject* LBField = PyList_New(0);
  PyObject* LTranslation = PyList_New(0);
  PyObject* LRotation = PyList_New(0);
  PyObject* LWidth = PyList_New(0);

  TVector3D BField(0, 0, 0);
  TVector3D Width (0, 0, 0);
  TVector3D Rotation(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  static char *kwlist[] = {"bfield", "width", "rotations", "translation", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|OO", kwlist,
                                                          &LBField,
                                                          &LWidth,
                                                          &LRotation,
                                                          &LTranslation)) {
    return NULL;
  }



  // Check BField
  if (PyList_Size(LBField) != 0) {
    BField = SRS_ListAsTVector3D(LBField);
  }


  // Check Width
  if (PyList_Size(LWidth) != 0) {
    Width = SRS_ListAsTVector3D(LWidth);
  }


  // Check Rotations
  if (PyList_Size(LRotation) != 0) {
    Rotation = SRS_ListAsTVector3D(LRotation);
  }


  // Check Translation
  if (PyList_Size(LTranslation) != 0) {
    Translation = SRS_ListAsTVector3D(LTranslation);
  }


  BField.RotateSelfXYZ(Rotation);
  Width.RotateSelfXYZ(Rotation);

  // Set the object variable
  self->obj->AddMagneticField((TBField*) new TBField3D_Uniform(BField, Width, Translation));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}





static PyObject* SRS_ClearMagneticFields (SRSObject* self)
{
  // Clear all magnetic fields in the SRS object

  self->obj->ClearMagneticFields();

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}










static PyObject* SRS_GetBField (SRSObject* self, PyObject* args)
{
  // Get the magnetic field at a point as a 3D list [Bx, By, Bz]

  // Python list object
  PyObject * List;

  // Grab the python list from the input
  if (! PyArg_ParseTuple( args, "O!", &PyList_Type, &List)) {
    return NULL;
  }

  // Has to have the correct number of arguments
  if (PyList_Size(List) != 3) {
    return NULL;
  }


  // Grab the values
  double const X = PyFloat_AsDouble(PyList_GetItem(List, 0));
  double const Y = PyFloat_AsDouble(PyList_GetItem(List, 1));
  double const Z = PyFloat_AsDouble(PyList_GetItem(List, 2));

  // Set the object variable
  TVector3D const B = self->obj->GetB(X, Y, Z);

  // Create a python list
  PyObject *PList = PyList_New(0);

  PyList_Append(PList, Py_BuildValue("f", B.GetX()));
  PyList_Append(PList, Py_BuildValue("f", B.GetY()));
  PyList_Append(PList, Py_BuildValue("f", B.GetZ()));

  // Return the python list
  return PList;
}








static PyObject* SRS_AddParticleBeam (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Add a particle beam to the experiment



  // Grab the values
  char* Type = "";
  char* Name = "";
  PyObject* LX0 = PyList_New(0);
  PyObject* LV0 = PyList_New(0);
  double Energy_GeV = 0;
  double T0 = 0;
  double Current = 0;
  double Weight = 1;
  double Mass = 0;
  double Charge = 0;
  PyObject* LRotation = PyList_New(0);
  PyObject* LTranslation = PyList_New(0);

  TVector3D X0;
  TVector3D V0;
  TVector3D Rotation(0, 0, 0);
  TVector3D Translation(0, 0, 0);



  static char *kwlist[] = {"type", "name", "x0", "v0", "energy_GeV", "t0", "current", "weight", "rotations", "translation", "mass", "charge", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ssOO|ddddOO", kwlist,
                                                               &Type,
                                                               &Name,
                                                               &LX0,
                                                               &LV0,
                                                               &Energy_GeV,
                                                               &T0,
                                                               &Current,
                                                               &Weight,
                                                               &LRotation,
                                                               &LTranslation,
                                                               &Mass,
                                                               &Charge)) {
    return NULL;
  }


  // Check that type and name exist
  if (Type == "" || Name == "") {
    throw;
  }

  // If this is a custom particle
  if (Type == "custom") {
    if (Mass == 0 || Charge == 0) {
      throw;
    }
    // UPDATE: for custom beams
    throw;
  }

  X0 = SRS_ListAsTVector3D(LX0);
  V0 = SRS_ListAsTVector3D(LV0);

  if (Energy_GeV == 0) {
    // UPDATE: get energy from V vector
  }

  // Check Rotations
  if (PyList_Size(LRotation) != 0) {
    Rotation = SRS_ListAsTVector3D(LRotation);
  }

  // Check Translation
  if (PyList_Size(LTranslation) != 0) {
    Translation = SRS_ListAsTVector3D(LTranslation);
  }

  // Add the particle beam
  self->obj->AddParticleBeam(Type, Name, X0, V0, Energy_GeV, T0, Current, Weight);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




















static PyObject* SRS_SetNewParticle (SRSObject* self)
{
  // Set a new particle within the SRS object

  self->obj->SetNewParticle();

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






static PyObject* SRS_ClearParticleBeams (SRSObject* self)
{
  // Clear the contents of the particle beam container in SRS

  self->obj->ClearParticleBeams();

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






static PyObject* SRS_CalculateTrajectory (SRSObject* self)
{
  // Get the CTStop variable from SRS

  self->obj->CalculateTrajectory();

  // Return the trajectory
  return SRS_GetTrajectory(self);
}






static PyObject* SRS_GetTrajectory (SRSObject* self)
{
  // Get the Trajectory as 2 3D lists [[x, y, z], [BetaX, BetaY, BetaZ]]

  // Create a python list
  PyObject *PList = PyList_New(0);

  // Grab trajectory
  TParticleTrajectoryPoints const& T = self->obj->GetTrajectory();

  // Number of points in trajectory calculation
  size_t NTPoints = T.GetNPoints();

  // Loop over all points in trajectory
  for (int iT = 0; iT != NTPoints; ++iT) {
    // Create a python list for X and Beta
    PyObject *PList2 = PyList_New(0);

    // Add position and Beta to list
    PyList_Append(PList2, SRS_TVector3DAsList(T.GetX(iT)));
    PyList_Append(PList2, SRS_TVector3DAsList(T.GetB(iT)));
    PyList_Append(PList, PList2);
  }

  // Return the python list
  return PList;
}






static PyObject* SRS_GetSpectrum (SRSObject* self)
{
  // Get the Trajectory as 2 3D lists [[x, y, z], [BetaX, BetaY, BetaZ]]

  // Create a python list
  PyObject *PList = PyList_New(0);

  // Grab trajectory
  TSpectrumContainer const& Spectrum = self->obj->GetSpectrum();

  // Number of points in trajectory calculation
  size_t NSPoints = Spectrum.GetNPoints();

  // Loop over all points in trajectory
  for (int iS = 0; iS != NSPoints; ++iS) {
    // Create a python list for X and Beta
    PyObject *PList2 = PyList_New(0);

    // Add position and Beta to list
    PyList_Append(PList2, Py_BuildValue("f", Spectrum.GetEnergy(iS)));
    PyList_Append(PList2, Py_BuildValue("f", Spectrum.GetFlux(iS)));
    PyList_Append(PList, PList2);
  }

  // Return the python list
  return PList;
}

























static PyObject* SRS_CalculateSpectrum (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Calculate the spectrum given an observation point, and energy range


  PyObject* LX0 = PyList_New(0);
  int NPoints = 0;
  PyObject* LEnergyRange_eV = PyList_New(0);
  PyObject* LPoints_eV = PyList_New(0);



  static char *kwlist[] = {"xyz", "npoints", "energy_range_eV", "points_eV", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOO", kwlist,
                                                          &LX0,
                                                          &NPoints,
                                                          &LEnergyRange_eV,
                                                          &LPoints_eV)) {
    return NULL;
  }


  // Add all values to a vector
  std::vector<double> VEnergy;
  for (int i = 0; i < PyList_Size(LPoints_eV); ++i) {
    VEnergy.push_back(PyFloat_AsDouble(PyList_GetItem(LPoints_eV, i)));
  }

  double EStart;
  double EStop;

  if (PyList_Size(LEnergyRange_eV) != 0) {
    if (PyList_Size(LEnergyRange_eV) == 2) {
      EStart = PyFloat_AsDouble(PyList_GetItem(LEnergyRange_eV, 0));
      EStop  = PyFloat_AsDouble(PyList_GetItem(LEnergyRange_eV, 1));
    } else {
      throw;
    }
  }




  // Observation point
  TVector3D const X0 = SRS_ListAsTVector3D(LX0);

  // Actually calculate the spectrum
  if (VEnergy.size() == 0) {
    self->obj->CalculateSpectrum(X0, EStart, EStop, NPoints);
  } else {
    self->obj->CalculateSpectrum(X0, VEnergy);
  }

  // Return the spectrum
  return SRS_GetSpectrum(self);
}









static PyObject* SRS_CalculateTotalPower (SRSObject* self)
{
  // Calculate the total power radiated by the current particle

  // Return the total power
  return Py_BuildValue("f", self->obj->CalculateTotalPower());
}

















static PyObject* SRS_CalculatePowerDensityRectangle (SRSObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  char*     SurfacePlane = "";
  size_t    NX1 = 0;
  size_t    NX2 = 0;
  double    Width_X1 = 0;
  double    Width_X2 = 0;
  PyObject* List_NPoints= PyList_New(0);
  PyObject* List_Width= PyList_New(0);
  PyObject* List_Translation = PyList_New(0);
  PyObject* List_Rotations = PyList_New(0);
  PyObject* List_X0X1X2 = PyList_New(0);
  int       NormalDirection = 0;
  int       Dim = 2;


  static char *kwlist[] = {"npoints", "plane", "normal", "dim", "width", "rotations", "translation", "x0x1x2", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|siiOOOO", kwlist,
                                                              &List_NPoints,
                                                              &SurfacePlane,
                                                              &NormalDirection,
                                                              &Dim,
                                                              &List_Width,
                                                              &List_Rotations,
                                                              &List_Translation,
                                                              &List_X0X1X2)) {
    return NULL;
  }


  // The rectangular surface object we'll use
  TSurfacePoints_Rectangle Surface;

  if (PyList_Size(List_NPoints) == 2) {
    // NPoints in [m]
    NX1 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 0));
    NX2 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 1));
  } else {
    throw;
  }



  if (NX1 <= 0 || NX2 <= 0) {
    std::cerr << "ERROR: NX1,2 < 1" << std::endl;
    throw;
  }

  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    Rotations = SRS_ListAsTVector3D(List_Rotations);
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    Translation = SRS_ListAsTVector3D(List_Translation);
  }

  if (PyList_Size(List_Width) == 2) {
    // Width in [m]
    Width_X1 = PyFloat_AsDouble(PyList_GetItem(List_Width, 0));
    Width_X2 = PyFloat_AsDouble(PyList_GetItem(List_Width, 1));
  }



  // If you are requesting a simple surface plane, check that you have widths
  if (SurfacePlane != "" && Width_X1 > 0 && Width_X2 > 0) {
    Surface.Init(SurfacePlane, NX1, NX2, Width_X1, Width_X2, Rotations, Translation, NormalDirection);
  }



  // If X0X1X2 defined
  std::vector<TVector3D> X0X1X2;

  if (PyList_Size(List_X0X1X2) != 0) {
    if (PyList_Size(List_X0X1X2) == 3) {
      for (int i = 0; i != 3; ++i) {
        PyObject* List_X = PyList_GetItem(List_X0X1X2, i);

        X0X1X2.push_back(SRS_ListAsTVector3D(List_X));
      }
    } else {
      throw;
    }

    // UPDATE: Check for orthogonality
    Surface.Init(NX1, NX2, X0X1X2[0], X0X1X2[1], X0X1X2[2], NormalDirection);
  }


  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;


  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;
  self->obj->CalculatePowerDensity(Surface, PowerDensityContainer, Dim, Directional);


  // Build the output list of: [[[x, y, z], PowerDensity], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = PowerDensityContainer.GetNPoints();

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = PowerDensityContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    PyList_Append(PList2, SRS_TVector3DAsList(P.GetX()));
    PyList_Append(PList2, Py_BuildValue("f", P.GetV()));
    PyList_Append(PList, PList2);

  }

  return PList;
}














static PyObject* SRS_CalculateFluxRectangle (SRSObject* self, PyObject* args)
{
  // Calculate the spectrum given an observation point, and energy range


  double    Energy;
  char*     SurfacePlane;
  int       NX1;
  int       NX2;
  double    Width_X1;
  double    Width_X2;
  PyObject* List_Observer;
  int       NormalDirection;


  // Grab the values
  if (! PyArg_ParseTuple(args, "dsdidiO!i", &Energy, &SurfacePlane, &Width_X1, &NX1, &Width_X2, &NX2, &PyList_Type, &List_Observer, &NormalDirection)) {
    return NULL;
  }

  // Has to have the correct number of arguments
  if (PyList_Size(List_Observer) != 3) {
    return NULL;
  }

  // Observation point
  TVector3D const ObservationPoint = SRS_ListAsTVector3D(List_Observer);

  // Container for Point plus scalar
  T3DScalarContainer FluxContainer;

  TSurfacePoints_RectangleSimple Surface(SurfacePlane, NX1, NX2, Width_X1, Width_X2, ObservationPoint, NormalDirection);

  // Actually calculate the spectrum
  self->obj->CalculateFlux(Surface, Energy, FluxContainer);


  // Build the output list of: [[[x, y, z], Flux], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = FluxContainer.GetNPoints();

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar F = FluxContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);

    // Add position and value to list
    PyList_Append(PList2, SRS_TVector3DAsList(F.GetX()));
    PyList_Append(PList2, Py_BuildValue("f", F.GetV()));
    PyList_Append(PList, PList2);

  }

  return PList;
}








































static PyGetSetDef SRS_getseters[] = {
  {"ctstart", (getter)SRS_GetCTStart, (setter)SRS_SetCTStart, "Start time for calculation in [m]", NULL},
  {"ctstop",  (getter)SRS_GetCTStop,  (setter)SRS_SetCTStop,  "Stop time for calculation in [m]", NULL},
  {"npoints_trajectory",  (getter)SRS_GetNPointsTrajectory,  (setter)SRS_SetNPointsTrajectory,  "Total number of points for trajectory", NULL},
  {NULL}  /* Sentinel */
};




static PyMethodDef SRS_methods[] = {
  // We must tell python about the function we allow access as well as give them nice
  // python names, and tell python the method of input parameters.

  {"pi",                       (PyCFunction) SRS_Pi,                      METH_NOARGS,  "do you want some pi?"},


  {"set_ctstart",                       (PyCFunction) SRS_SetCTStart,                      METH_O,       
     "set the start time in [m]\n\n\\
     :param name: the input name and a much longer description of what this function does really\n\\
     :type name: str\n\\
     :param last: the input last name\n\\
     :type last: str\n\\
     :returns: str -- a string hello\n\\
     :raises: AttributeError\n"
  },
  //{"set_ctstop2",                       (PyCFunction) SRS_SetCTStop2,                      METH_VARARGS | METH_KEYWORDS,       "set the stop time in [m]"},
  {"get_ctstart",                       (PyCFunction) SRS_GetCTStart,                      METH_NOARGS,  "get the start time in [m]"},
  {"set_ctstop",                        (PyCFunction) SRS_SetCTStop,                       METH_O,       "set the stop time in [m]"},
  {"get_ctstop",                        (PyCFunction) SRS_GetCTStop,                       METH_NOARGS,  "get the stop time in [m]"},
  {"set_npoints_trajectory",            (PyCFunction) SRS_SetNPointsTrajectory,            METH_O,       "set the total number of points for the trajectory"},
  {"get_npoints_trajectory",            (PyCFunction) SRS_GetNPointsTrajectory,            METH_NOARGS,  "get the total number of points for the trajectory"},
  {"set_ctstartstop",                   (PyCFunction) SRS_SetCTStartStop,                  METH_VARARGS, "set the start and stop time in [m]"},
                                                                                          
  {"add_magnetic_field",                (PyCFunction) SRS_AddMagneticField,                METH_VARARGS | METH_KEYWORDS, "add a magnetic field from a file"},
  {"add_magnetic_field_function",       (PyCFunction) SRS_AddMagneticFieldFunction,        METH_VARARGS, "add a magnetic field in form of python function"},
  {"add_magnetic_field_gaussian",       (PyCFunction) SRS_AddMagneticFieldGaussian,        METH_VARARGS | METH_KEYWORDS, "add a magnetic field in form of 3D gaussian"},
  {"clear_magnetic_fields",             (PyCFunction) SRS_ClearMagneticFields,             METH_NOARGS,  "clear all internal magnetic fields"},
  {"add_magnetic_field_uniform",        (PyCFunction) SRS_AddMagneticFieldUniform,         METH_VARARGS | METH_KEYWORDS, "add a uniform magnetic field in 3D"},
  {"get_bfield",                        (PyCFunction) SRS_GetBField,                       METH_VARARGS, "get the magnetic field at a given position in space (and someday time?)"},
                                                                                          
                                                                                          
  {"add_particle_beam",                 (PyCFunction) SRS_AddParticleBeam,                 METH_VARARGS | METH_KEYWORDS, "add a particle beam"},
  {"clear_particle_beams",              (PyCFunction) SRS_ClearParticleBeams,              METH_NOARGS,  "Clear all existing particle beams from SRS"},
                                                                                          
  {"set_new_particle",                  (PyCFunction) SRS_SetNewParticle,                  METH_NOARGS,  "Set the internal particle to a new random particle"},
                                                                                          
  {"calculate_trajectory",              (PyCFunction) SRS_CalculateTrajectory,             METH_NOARGS,  "Calclate the trajectory for the current particle"},
  {"get_trajectory",                    (PyCFunction) SRS_GetTrajectory,                   METH_NOARGS,  "Get the trajectory for the current particle as 2 3D lists [[x, y, z], [BetaX, BetaY, BetaZ]]"},

  {"calculate_spectrum",                (PyCFunction) SRS_CalculateSpectrum,               METH_VARARGS | METH_KEYWORDS, "calculate the spectrum at an observation point"},
  {"get_spectrum",                      (PyCFunction) SRS_GetSpectrum,                     METH_NOARGS,  "Get the spectrum for the current particle as list of 2D [[energy, flux], [...]...]"},

  {"calculate_total_power",             (PyCFunction) SRS_CalculateTotalPower,             METH_NOARGS,  "calculate total power radiated"},
  {"calculate_power_density_rectangle", (PyCFunction) SRS_CalculatePowerDensityRectangle,  METH_VARARGS | METH_KEYWORDS, "calculate the power density given a surface"},
  {"calculate_flux_rectangle",          (PyCFunction) SRS_CalculateFluxRectangle,          METH_VARARGS, "calculate the flux given a surface"},

  {NULL}  /* Sentinel */
};




static PyTypeObject SRSType = {
  // The python object.  Fully defined elsewhere.  only put here what you need,
  // otherwise default values

  PyObject_HEAD_INIT(NULL)
  0,                                        /* ob_size */
  "mod.SRS",                                 /* tp_name */
  sizeof(SRSObject),                         /* tp_basicsize */
  0,                                        /* tp_itemsize */
  (destructor) SRS_dealloc,                 /* tp_dealloc */
  0,                                        /* tp_print */
  0,                                        /* tp_getattr */
  0,                                        /* tp_setattr */
  0,                                        /* tp_compare */
  0,                                        /* tp_repr */
  0,                                        /* tp_as_number */
  0,                                        /* tp_as_sequence */
  0,                                        /* tp_as_mapping */
  0,                                        /* tp_hash */
  0,                                        /* tp_call */
  0,                                        /* tp_str */
  0,                                        /* tp_getattro */
  0,                                        /* tp_setattro */
  0,                                        /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
  "SRS class",                              /* tp_doc */
  0,                                        /* tp_traverse */
  0,                                        /* tp_clear */
  0,                                        /* tp_richcompare */
  0,                                        /* tp_weaklistoffset */
  0,                                        /* tp_iter */
  0,                                        /* tp_iternext */
  SRS_methods,                             /* tp_methods */
  0,                                        /* tp_members */
  SRS_getseters,                           /* tp_getset */
  0,                                        /* tp_base */
  0,                                        /* tp_dict */
  0,                                        /* tp_descr_get */
  0,                                        /* tp_descr_set */
  0,                                        /* tp_dictoffset */
  0,                                        /* tp_init */
  0,                                        /* tp_alloc */
  SRS_new,                                  /* tp_new */
};




static PyMethodDef module_methods[] = {
  // I do not need
  {NULL}  /* Sentinel */
};




PyMODINIT_FUNC initSRS ()
{
  // Initialization of the module.


  PyObject* m;

  if (PyType_Ready(&SRSType) < 0)
    return;

  m = Py_InitModule3("SRS", module_methods, "SRS Extension module for the SRS Python API.");
  if (m == NULL)
    return;

  Py_INCREF(&SRSType);
  PyModule_AddObject(m, "SRS", (PyObject *)&SRSType);
}


