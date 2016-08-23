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

#include "TSurfacePoints_Rectangle.h"
#include "TSurfacePoints_3D.h"
#include "T3DScalarContainer.h"
#include "TBFieldPythonFunction.h"
#include "TBField3D_Gaussian.h"
#include "TBField3D_Uniform.h"
#include "TBField3D_IdealUndulator.h"

#include <iostream>
#include <vector>
#include <algorithm>





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







static TVector2D SRS_ListAsTVector2D (PyObject* List)
{
  TVector2D V;
  if (PyList_Size(List) == 2) {
    Py_INCREF(List);
    V.SetXY(PyFloat_AsDouble(PyList_GetItem(List, 0)),
             PyFloat_AsDouble(PyList_GetItem(List, 1)));
    Py_DECREF(List);
  } else {
    throw std::length_error("number of elements not 2");
  }

  // Return the python list
  return V;
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
    throw std::length_error("number of elements not 3");
  }

  // Return the python list
  return V;
}







static PyObject* SRS_TVector3DAsList (TVector3D const& V)
{
  // Turn a TVector3D into a list (like a vector)

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
  // Return the internal SRS number constant pi
  return Py_BuildValue("d", TSRS::Pi());
}



static PyObject* SRS_Qe (SRSObject* self, PyObject* arg)
{
  // Return the internal SRS number for elementary charge
  return Py_BuildValue("d", TSRS::Qe());
}




static PyObject* SRS_Me (SRSObject* self, PyObject* arg)
{
  // Return the internal SRS number for mass of the electron [kg]
  return Py_BuildValue("d", TSRS::Me());
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
  if (!PyArg_ParseTuple(args, "dd", &Start, &Stop)) {
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
  // Add a magnetic field from a file.
  // UPDATE: add binary file reading


  // Grab the values
  char const* FileName = "";
  char const* FileFormat = "";
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Scaling     = PyList_New(0);

  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  std::vector<double> Scaling;


  // Input variables and parsing
  static char *kwlist[] = {"ifile", "iformat", "rotations", "translation", "scale", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ss|OOO", kwlist,
                                                           &FileName,
                                                           &FileFormat,
                                                           &List_Rotations,
                                                           &List_Translation,
                                                           &List_Scaling)) {
    return NULL;
  }

  // Check that filename and format exist
  if (FileName == "" || FileFormat == "") {
    PyErr_SetString(PyExc_ValueError, "'ifile' or 'iformat' is blank");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Get any scaling factors
  // UPDATE: Check against fileformat number of strings
  for (int i = 0; i < PyList_Size(List_Scaling); ++i) {
    Scaling.push_back(PyFloat_AsDouble(PyList_GetItem(List_Scaling, i)));
  }

  // Add the magnetic field to the SRS object
  try {
    self->obj->AddMagneticField(FileName, FileFormat, Rotations, Translation, Scaling);
  } catch (...) {
    PyErr_SetString(PyExc_ValueError, "Could not import magnetic field.  Check 'ifile' and 'iformat' are correct");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRS_AddMagneticFieldFunction (SRSObject* self, PyObject* args)
{
  // Add a python function as a magnetic field object

  // Grab the function
  PyObject* Function;
  if (! PyArg_ParseTuple(args, "O:set_callback", &Function)) {
    return NULL;
  }

  // Increment ref to function for python
  Py_INCREF(Function);

  // Add the function as a field to the SRS object
  try {
    self->obj->AddMagneticField( (TBField*) new TBFieldPythonFunction(Function));
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  // Decrement reference to function for python
  Py_DECREF(Function);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRS_AddMagneticFieldGaussian (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Add a magnetic field that is a gaussian

  // Lists and variables
  PyObject* List_BField       = PyList_New(0);
  PyObject* List_Translation  = PyList_New(0);
  PyObject* List_Rotations    = PyList_New(0);
  PyObject* List_Sigma        = PyList_New(0);

  TVector3D BField(0, 0, 0);
  TVector3D Sigma(0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Input variables and parsing
  static char *kwlist[] = {"bfield", "sigma", "rotations", "translation", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OO|OO", kwlist,
                                                          &List_BField,
                                                          &List_Sigma,
                                                          &List_Rotations,
                                                          &List_Translation)) {
    return NULL;
  }



  // Check BField
  try {
    BField = SRS_ListAsTVector3D(List_BField);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield'");
    return NULL;
  }

  // Check Width
  try {
    Sigma = SRS_ListAsTVector3D(List_Sigma);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'sigma'");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Rotate field and sigma
  // UPDATE: check this
  BField.RotateSelfXYZ(Rotations);
  Sigma.RotateSelfXYZ(Rotations);

  // Add field
  self->obj->AddMagneticField( (TBField*) new TBField3D_Gaussian(BField, Translation, Sigma));

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}








static PyObject* SRS_AddMagneticFieldUniform (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Add a uniform field with a given width in a given direction, or for all space

  // Lists and vectors
  PyObject* List_BField      = PyList_New(0);
  PyObject* List_Translation = PyList_New(0);
  PyObject* List_Rotations   = PyList_New(0);
  PyObject* List_Width       = PyList_New(0);

  TVector3D BField(0, 0, 0);
  TVector3D Width (0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Input variables and parsing
  static char *kwlist[] = {"bfield", "width", "rotations", "translation", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|OOO", kwlist,
                                                          &List_BField,
                                                          &List_Width,
                                                          &List_Rotations,
                                                          &List_Translation)) {
    return NULL;
  }


  // Check BField
  try {
    BField = SRS_ListAsTVector3D(List_BField);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield'");
    return NULL;
  }

  // Check Width
  if (PyList_Size(List_Width) != 0) {
    try {
      Width = SRS_ListAsTVector3D(List_Width);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'width'");
      return NULL;
    }
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  // Rotate field and widths
  BField.RotateSelfXYZ(Rotations);
  Width.RotateSelfXYZ(Rotations);

  // Add the field
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










static PyObject* SRS_AddMagneticFieldIdealUndulator (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Set the start and stop times for SRS in [m]

  // Lists and variables
  PyObject* List_BField      = PyList_New(0);
  PyObject* List_Period      = PyList_New(0);
  PyObject* List_Rotations   = PyList_New(0);
  PyObject* List_Translation = PyList_New(0);
  int       NPeriods = 0;
  double    Phase = 0;

  TVector3D BField(0, 0, 0);
  TVector3D Period(0, 0, 0);
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Input variables and parsing
  static char *kwlist[] = {"bfield", "period", "nperiods", "phase", "rotations", "translation", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOi|dOO", kwlist,
                                                            &List_BField,
                                                            &List_Period,
                                                            &NPeriods,
                                                            &Phase,
                                                            &List_Rotations,
                                                            &List_Translation)) {
    return NULL;
  }


  // Check BField
  try {
    BField = SRS_ListAsTVector3D(List_BField);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'bfield'");
    return NULL;
  }

  // Check Period
  try {
    Period = SRS_ListAsTVector3D(List_Period);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'period'");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }

  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Rotate field and sigma
  // UPDATE: check this
  BField.RotateSelfXYZ(Rotations);
  Period.RotateSelfXYZ(Rotations);


  // Add field
  self->obj->AddMagneticField( (TBField*) new TBField3D_IdealUndulator(BField, Period, NPeriods, Translation, Phase));

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
  TVector3D X;
  try {
    X = SRS_ListAsTVector3D(List);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
    return NULL;
  }

  // Set the object variable
  TVector3D const B = self->obj->GetB(X);

  // Create a python list
  PyObject *PList = SRS_TVector3DAsList(B);

  // Return the python list
  return PList;
}



static PyObject* SRS_SetParticleBeam (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Clear all particle beams, add this beam, and set a new particle

  self->obj->ClearParticleBeams();

  SRS_AddParticleBeam(self, args, keywds);

  self->obj->SetNewParticle("", "ideal");

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}





static PyObject* SRS_AddParticleBeam (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Add a particle beam to the experiment

  // Lists and variables some with initial values
  char const* Type                       = "";
  char const* Name                       = "";
  double      Energy_GeV                 = 0;
  double      Sigma_Energy_GeV           = 0;
  double      T0                         = 0;
  double      Current                    = 0;
  double      Weight                     = 1;
  double      Mass                       = 0;
  double      Charge                     = 0;
  PyObject*   List_X0                    = PyList_New(0);
  PyObject*   List_Direction             = PyList_New(0);
  PyObject*   List_Rotations             = PyList_New(0);
  PyObject*   List_Translation           = PyList_New(0);
  PyObject*   List_Horizontal_Direction  = PyList_New(0);
  PyObject*   List_Beta                  = PyList_New(0);
  PyObject*   List_Emittance             = PyList_New(0);
  PyObject*   List_Lattice_Center        = PyList_New(0);

  TVector3D X0(0, 0, 0);
  TVector3D Direction;
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);
  TVector3D Horizontal_Direction;
  TVector2D Beta(0, 0);
  TVector2D Emittance(0, 0);
  TVector3D Lattice_Center(0, 0, 0);


  // Input variables and parsing
  static char *kwlist[] = {"type", "name", "energy_GeV", "direction", "x0", "sigma_energy_GeV", "t0", "current", "weight", "rotations", "translation", "horizontal_direction", "beta", "emittance", "lattice_center", "mass", "charge", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "ssdO|OddddOOOOOOdd", kwlist,
                                                                       &Type,
                                                                       &Name,
                                                                       &Energy_GeV,
                                                                       &List_Direction,
                                                                       &List_X0,
                                                                       &Sigma_Energy_GeV,
                                                                       &T0,
                                                                       &Current,
                                                                       &Weight,
                                                                       &List_Rotations,
                                                                       &List_Translation,
                                                                       &List_Horizontal_Direction,
                                                                       &List_Beta,
                                                                       &List_Emittance,
                                                                       &List_Lattice_Center,
                                                                       &Mass,
                                                                       &Charge)) {
    return NULL;
  }


  // Check that type and name exist
  if (Type == "" || Name == "") {
    PyErr_SetString(PyExc_ValueError, "'type' or 'name' is blank");
    return NULL;
  }

  // If this is a custom particle
  if (PyList_Size(List_X0) != 0) {
    try {
      X0 = SRS_ListAsTVector3D(List_X0);
      Direction = SRS_ListAsTVector3D(List_Direction);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0' or 'direction'");
      return NULL;
    }
  }

  if (Sigma_Energy_GeV == 0) {
    // Do nothing.  zero energy diff is alright
  } else if (Sigma_Energy_GeV < 0) {
    PyErr_SetString(PyExc_ValueError, "'sigma_energy_GeV' cannot be less than zero");
    return NULL;
  }

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }



  // Check for Horizontal_Direction in the input
  if (PyList_Size(List_Horizontal_Direction) != 0) {
    try {
      Horizontal_Direction = SRS_ListAsTVector3D(List_Horizontal_Direction);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'horizontal_direction'");
      return NULL;
    }
  } else {
    Horizontal_Direction = Direction.Orthogonal().UnitVector();
  }
  Horizontal_Direction = Horizontal_Direction.UnitVector();


  // Check for Beta in the input
  if (PyList_Size(List_Beta) != 0) {
    try {
      Beta = SRS_ListAsTVector2D(List_Beta);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'beta'");
      return NULL;
    }
  }


  // Check for Emittance in the input
  if (PyList_Size(List_Emittance) != 0) {
    try {
      Emittance = SRS_ListAsTVector2D(List_Emittance);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'emittance'");
      return NULL;
    }
  }


  // Check for Lattice Center in the input
  if (PyList_Size(List_Lattice_Center) != 0) {
    try {
      Lattice_Center = SRS_ListAsTVector3D(List_Lattice_Center);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'lattice_center'");
      return NULL;
    }
  }


  // If emittance and beta are defined get RMS values
  TVector2D SigmaU(0, 0);
  TVector2D SigmaUPrime(0, 0);
  if (PyList_Size(List_Emittance) != 0 && PyList_Size(List_Beta) != 0) {
    SigmaU[0]      = sqrt(Emittance[0] * Beta[0]);
    SigmaU[1]      = sqrt(Emittance[1] * Beta[1]);
    SigmaUPrime[0] = sqrt(Emittance[0] / Beta[0]);
    SigmaUPrime[1] = sqrt(Emittance[1] / Beta[1]);
  }




  // Rotate beam parameters
  X0.RotateSelfXYZ(Rotations);
  X0 += Translation;

  // UPDATE: An idea for later, use new variable "velocity"
  //if (Energy_GeV == 0) {
  //  if (Direction.Mag() >= 1) {
  //    throw;
  //  }
  //  Energy_GeV = sqrt(1.0 / (1.0 - Direction.Mag2())) * Mass;
  //}


  // Add the particle beam
  try {
    if (std::string(Type) == "custom") {
      if (Mass == 0 || Charge == 0) {
        PyErr_SetString(PyExc_ValueError, "'mass' or 'charge' is zero");
        return NULL;
      }
      // UPDATE: for custom beams
      self->obj->AddParticleBeam(Type, Name, X0, Direction, Energy_GeV, T0, Current, Weight, Charge, Mass);
    } else {
      self->obj->AddParticleBeam(Type, Name, X0, Direction, Energy_GeV, T0, Current, Weight);
    }
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, "invalid argument in adding particle beam.  possibly 'name' already exists");
    return NULL;
  }

  // UPDATE: Change me
  self->obj->GetParticleBeam(Name).SetSigma(Horizontal_Direction, SigmaU, SigmaUPrime, Lattice_Center, Sigma_Energy_GeV);

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}




















static PyObject* SRS_SetNewParticle (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Set a new particle within the SRS object

  char const* Beam_IN = "";
  char const* Particle_IN = "";

  static char *kwlist[] = {"beam", "particle", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ss", kwlist, &Beam_IN, &Particle_IN)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }

  std::string Beam = Beam_IN;
  std::string Particle = Particle_IN;

  std::transform(Beam.begin(), Beam.end(), Beam.begin(), ::tolower);
  std::transform(Particle.begin(), Particle.end(), Particle.begin(), ::tolower);

  if (Particle != "" && !(Particle == "ideal" || Particle == "random")) {
    PyErr_SetString(PyExc_ValueError, "'particle' must be 'ideal' or 'random'");
    return NULL;
  }

  try {
    if ((Beam) == "" && (Particle) == "") {
      self->obj->SetNewParticle();
    } else if (Beam != "" && Particle != "") {
      self->obj->SetNewParticle(Beam, Particle);
    } else if (Beam == "" && Particle != "") {
      self->obj->SetNewParticle(Beam, Particle);
    } else if (Beam != "" && Particle == "") {
      self->obj->SetNewParticle(Beam, Particle);
    }
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, "'beam' name not found");
    return NULL;
  }

  // Must return python object None in a special way
  Py_INCREF(Py_None);
  return Py_None;
}






static PyObject* SRS_GetParticleX0 (SRSObject* self)
{
  // Get the particle position at particle t0
  return SRS_TVector3DAsList( self->obj->GetCurrentParticle().GetX0() );
}




static PyObject* SRS_GetParticleBeta0 (SRSObject* self)
{
  // Get the particle beta at particle t0
  return SRS_TVector3DAsList( self->obj->GetCurrentParticle().GetB0() );
}




static PyObject* SRS_GetParticleE0 (SRSObject* self)
{
  // Get the particle beta at particle t0
  return Py_BuildValue("f", (self->obj->GetCurrentParticle().GetE0()));
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

  try {
    self->obj->CalculateTrajectory();
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

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












static PyObject* SRS_GetSpectrumAsList (SRSObject* self, TSpectrumContainer const& Spectrum)
{
  // Get the Trajectory as 2 3D lists [[x, y, z], [BetaX, BetaY, BetaZ]]

  // Create a python list
  PyObject *PList = PyList_New(0);

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


  PyObject* List_Obs = PyList_New(0);
  int NPoints = 0;
  PyObject* List_EnergyRange_eV = PyList_New(0);
  PyObject* List_Points_eV      = PyList_New(0);
  int NParticles = 0;
  char* OFile = "";



  static char *kwlist[] = {"obs", "npoints", "energy_range_eV", "points_eV", "nparticles", "ofile", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOOis", kwlist,
                                                            &List_Obs,
                                                            &NPoints,
                                                            &List_EnergyRange_eV,
                                                            &List_Points_eV,
                                                            &NParticles,
                                                            &OFile)) {
    return NULL;
  }


  // Check if a beam is at least defined
  //if (self->obj->GetNParticleBeams() < 1) {
  //  PyErr_SetString(PyExc_ValueError, "No particle beam defined");
  //  return NULL;
  //}


  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }


  // Add all values to a vector
  std::vector<double> VPoints_eV;
  for (int i = 0; i < PyList_Size(List_Points_eV); ++i) {
    VPoints_eV.push_back(PyFloat_AsDouble(PyList_GetItem(List_Points_eV, i)));
  }

  double EStart;
  double EStop;

  if (PyList_Size(List_EnergyRange_eV) != 0) {
    if (PyList_Size(List_EnergyRange_eV) == 2) {
      EStart = PyFloat_AsDouble(PyList_GetItem(List_EnergyRange_eV, 0));
      EStop  = PyFloat_AsDouble(PyList_GetItem(List_EnergyRange_eV, 1));
    } else {
      PyErr_SetString(PyExc_ValueError, "'energy_range_eV' must be a list of length 2");
      return NULL;
    }
  }




  // Observation point
  TVector3D Obs(0, 0, 0);
  try {
    Obs = SRS_ListAsTVector3D(List_Obs);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'obs'");
    return NULL;
  }



  TSpectrumContainer SpectrumContainer;

  if (VPoints_eV.size() == 0) {
    SpectrumContainer.Init(NPoints, EStart, EStop);
  } else {
    SpectrumContainer.Init(VPoints_eV);
  }

  // Actually calculate the spectrum
  try {
    if (NParticles == 0) {
      self->obj->CalculateSpectrum(Obs, SpectrumContainer);
    } else {
      double const Weight = 1.0 / (double) NParticles;
      for (int i = 0; i != NParticles; ++i) {
        self->obj->SetNewParticle();
        self->obj->CalculateSpectrum(Obs, SpectrumContainer, Weight);
      }
    }
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }


  if (std::string(OFile) != "") {
    SpectrumContainer.SaveToFile(OFile);
  }

  // Return the spectrum
  return SRS_GetSpectrumAsList(self, SpectrumContainer);
}






static PyObject* SRS_CalculateSpectrumGPU (SRSObject* self, PyObject* args, PyObject* keywds)
{
  // Calculate the spectrum given an observation point, and energy range


  PyObject* List_Obs = PyList_New(0);
  int NPoints = 0;
  PyObject* List_EnergyRange_eV = PyList_New(0);
  PyObject* List_Points_eV      = PyList_New(0);
  int NParticles = 0;
  char* OFile = "";



  static char *kwlist[] = {"obs", "npoints", "energy_range_eV", "points_eV", "nparticles", "ofile", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOOis", kwlist,
                                                            &List_Obs,
                                                            &NPoints,
                                                            &List_EnergyRange_eV,
                                                            &List_Points_eV,
                                                            &NParticles,
                                                            &OFile)) {
    return NULL;
  }


  // Check if a beam is at least defined
  //if (self->obj->GetNParticleBeams() < 1) {
  //  PyErr_SetString(PyExc_ValueError, "No particle beam defined");
  //  return NULL;
  //}


  // Check number of particles
  if (NParticles < 0) {
    PyErr_SetString(PyExc_ValueError, "'nparticles' must be >= 1 (sort of)");
    return NULL;
  }


  // Add all values to a vector
  std::vector<double> VPoints_eV;
  for (int i = 0; i < PyList_Size(List_Points_eV); ++i) {
    VPoints_eV.push_back(PyFloat_AsDouble(PyList_GetItem(List_Points_eV, i)));
  }

  double EStart;
  double EStop;

  if (PyList_Size(List_EnergyRange_eV) != 0) {
    if (PyList_Size(List_EnergyRange_eV) == 2) {
      EStart = PyFloat_AsDouble(PyList_GetItem(List_EnergyRange_eV, 0));
      EStop  = PyFloat_AsDouble(PyList_GetItem(List_EnergyRange_eV, 1));
    } else {
      PyErr_SetString(PyExc_ValueError, "'energy_range_eV' must be a list of length 2");
      return NULL;
    }
  }




  // Observation point
  TVector3D Obs(0, 0, 0);
  try {
    Obs = SRS_ListAsTVector3D(List_Obs);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'obs'");
    return NULL;
  }



  TSpectrumContainer SpectrumContainer;

  if (VPoints_eV.size() == 0) {
    SpectrumContainer.Init(NPoints, EStart, EStop);
  } else {
    SpectrumContainer.Init(VPoints_eV);
  }

  // Actually calculate the spectrum
  try {
    if (NParticles == 0) {
      self->obj->CalculateSpectrumGPU(Obs, SpectrumContainer);
    } else {
      double const Weight = 1.0 / (double) NParticles;
      for (int i = 0; i != NParticles; ++i) {
        self->obj->SetNewParticle();
        self->obj->CalculateSpectrumGPU(Obs, SpectrumContainer, Weight);
      }
    }
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }


  if (std::string(OFile) != "") {
    SpectrumContainer.SaveToFile(OFile);
  }

  // Return the spectrum
  return SRS_GetSpectrumAsList(self, SpectrumContainer);
}











static PyObject* SRS_CalculateTotalPower (SRSObject* self)
{
  // Calculate the total power radiated by the current particle

  double Power = 0;

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }

  // Return the total power
  // UPDATE: This does not fail when no beam defined
  try {
    Power = self->obj->CalculateTotalPower();
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  } catch (std::out_of_range e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }

  return Py_BuildValue("f", Power);
}








static PyObject* SRS_CalculatePowerDensity (SRSObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_Points      = PyList_New(0);
  int         NormalDirection = 0;
  int const   Dim = 3;
  char const* OutFileName = "";


  static char *kwlist[] = {"points", "normal", "rotations", "translation", "ofile", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOOs", kwlist,
                                                           &List_Points,
                                                           &NormalDirection,
                                                           &List_Rotations,
                                                           &List_Translation,
                                                           &OutFileName)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 & Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);


  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }


  // Look for arbitrary shape 3D points
  TSurfacePoints_3D Surface;
  for (int i = 0; i < PyList_Size(List_Points); ++i) {
    PyObject* LXN = PyList_GetItem(List_Points, i);
    TVector3D X;
    TVector3D N;
    if (PyList_Size(LXN) == 2) {

      try {
        X = SRS_ListAsTVector3D(PyList_GetItem(LXN, 0));
        N = SRS_ListAsTVector3D(PyList_GetItem(LXN, 1));
      } catch (std::length_error e) {
        PyErr_SetString(PyExc_ValueError, "Incorrect format in 'points': Point or Normal does not have 3 elements");
        return NULL;
      }

      // Rotate point and normal
      X.RotateSelfXYZ(Rotations);
      N.RotateSelfXYZ(Rotations);

      // Translate point, normal does not get translated
      X += Translation;

      Surface.AddPoint(X, N);
    } else {
      // input format error
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'points'");
      return NULL;
    }
  }



  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;


  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;

  self->obj->CalculatePowerDensity(Surface, PowerDensityContainer, Dim, Directional, OutFileName);

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









static PyObject* SRS_CalculatePowerDensityRectangle (SRSObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  char const* SurfacePlane = "";
  size_t      NX1 = 0;
  size_t      NX2 = 0;
  double      Width_X1 = 0;
  double      Width_X2 = 0;
  PyObject*   List_NPoints     = PyList_New(0);
  PyObject*   List_Width       = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_X0X1X2      = PyList_New(0);
  int         NormalDirection = 0;
  int         Dim = 2;
  char*       OutFileName = "";


  static char *kwlist[] = {"npoints", "plane", "width", "x0x1x2", "rotations", "translation", "ofile", "normal", "dim",  NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|sOOOOsii", kwlist,
                                                               &List_NPoints,
                                                               &SurfacePlane,
                                                               &List_Width,
                                                               &List_X0X1X2,
                                                               &List_Rotations,
                                                               &List_Translation,
                                                               &OutFileName,
                                                               &NormalDirection,
                                                               &Dim)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 & Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // The rectangular surface object we'll use
  TSurfacePoints_Rectangle Surface;

  if (PyList_Size(List_NPoints) == 2) {
    // NPoints in [m]
    NX1 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 0));
    NX2 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 1));
  } else {
    PyErr_SetString(PyExc_ValueError, "'npoints' must be [int, int]");
    return NULL;
  }



  if (NX1 <= 0 || NX2 <= 0) {
    PyErr_SetString(PyExc_ValueError, "an entry in 'npoints' is <= 0");
    return NULL;
  }

  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  if (PyList_Size(List_Width) == 2) {
    // Width in [m]
    Width_X1 = PyFloat_AsDouble(PyList_GetItem(List_Width, 0));
    Width_X2 = PyFloat_AsDouble(PyList_GetItem(List_Width, 1));
  }



  // If you are requesting a simple surface plane, check that you have widths
  if (SurfacePlane != "" && Width_X1 > 0 && Width_X2 > 0) {
    try {
      Surface.Init(SurfacePlane, NX1, NX2, Width_X1, Width_X2, Rotations, Translation, NormalDirection);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  }



  // If X0X1X2 defined
  std::vector<TVector3D> X0X1X2;

  if (PyList_Size(List_X0X1X2) != 0) {
    if (PyList_Size(List_X0X1X2) == 3) {
      for (int i = 0; i != 3; ++i) {
        PyObject* List_X = PyList_GetItem(List_X0X1X2, i);

        try {
          X0X1X2.push_back(SRS_ListAsTVector3D(List_X));
        } catch (std::length_error e) {
          PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0x1x2'");
          return NULL;
        }
      }
    } else {
      PyErr_SetString(PyExc_ValueError, "'x0x1x2' must have 3 XYZ points defined correctly");
      return NULL;
    }

    for (std::vector<TVector3D>::iterator it = X0X1X2.begin(); it != X0X1X2.end(); ++it) {
      it->RotateSelfXYZ(Rotations);
      *it += Translation;
    }

    // UPDATE: Check for orthogonality
    Surface.Init(NX1, NX2, X0X1X2[0], X0X1X2[1], X0X1X2[2], NormalDirection);
  }


  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;


  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;
  self->obj->CalculatePowerDensity(Surface, PowerDensityContainer, Dim, Directional, OutFileName);


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




static PyObject* SRS_CalculatePowerDensityRectangleGPU (SRSObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  char const* SurfacePlane = "";
  size_t      NX1 = 0;
  size_t      NX2 = 0;
  double      Width_X1 = 0;
  double      Width_X2 = 0;
  PyObject*   List_NPoints     = PyList_New(0);
  PyObject*   List_Width       = PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations   = PyList_New(0);
  PyObject*   List_X0X1X2      = PyList_New(0);
  int         NormalDirection = 0;
  int         Dim = 2;
  char*       OutFileName = "";


  static char *kwlist[] = {"npoints", "plane", "width", "x0x1x2", "rotations", "translation", "ofile", "normal", "dim",  NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|sOOOOsii", kwlist,
                                                               &List_NPoints,
                                                               &SurfacePlane,
                                                               &List_Width,
                                                               &List_X0X1X2,
                                                               &List_Rotations,
                                                               &List_Translation,
                                                               &OutFileName,
                                                               &NormalDirection,
                                                               &Dim)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 & Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // The rectangular surface object we'll use
  TSurfacePoints_Rectangle Surface;

  if (PyList_Size(List_NPoints) == 2) {
    // NPoints in [m]
    NX1 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 0));
    NX2 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 1));
  } else {
    PyErr_SetString(PyExc_ValueError, "'npoints' must be [int, int]");
    return NULL;
  }



  if (NX1 <= 0 || NX2 <= 0) {
    PyErr_SetString(PyExc_ValueError, "an entry in 'npoints' is <= 0");
    return NULL;
  }

  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  if (PyList_Size(List_Width) == 2) {
    // Width in [m]
    Width_X1 = PyFloat_AsDouble(PyList_GetItem(List_Width, 0));
    Width_X2 = PyFloat_AsDouble(PyList_GetItem(List_Width, 1));
  }



  // If you are requesting a simple surface plane, check that you have widths
  if (SurfacePlane != "" && Width_X1 > 0 && Width_X2 > 0) {
    try {
      Surface.Init(SurfacePlane, NX1, NX2, Width_X1, Width_X2, Rotations, Translation, NormalDirection);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  }



  // If X0X1X2 defined
  std::vector<TVector3D> X0X1X2;

  if (PyList_Size(List_X0X1X2) != 0) {
    if (PyList_Size(List_X0X1X2) == 3) {
      for (int i = 0; i != 3; ++i) {
        PyObject* List_X = PyList_GetItem(List_X0X1X2, i);

        try {
          X0X1X2.push_back(SRS_ListAsTVector3D(List_X));
        } catch (std::length_error e) {
          PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0x1x2'");
          return NULL;
        }
      }
    } else {
      PyErr_SetString(PyExc_ValueError, "'x0x1x2' must have 3 XYZ points defined correctly");
      return NULL;
    }

    for (std::vector<TVector3D>::iterator it = X0X1X2.begin(); it != X0X1X2.end(); ++it) {
      it->RotateSelfXYZ(Rotations);
      *it += Translation;
    }

    // UPDATE: Check for orthogonality
    Surface.Init(NX1, NX2, X0X1X2[0], X0X1X2[1], X0X1X2[2], NormalDirection);
  }


  // Container for Point plus scalar
  T3DScalarContainer PowerDensityContainer;


  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;
  try {
    self->obj->CalculatePowerDensityGPU(Surface, PowerDensityContainer, Dim, Directional, OutFileName);
  } catch (std::invalid_argument e) {
    PyErr_SetString(PyExc_ValueError, e.what());
    return NULL;
  }


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

















static PyObject* SRS_CalculateFluxRectangle (SRSObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the spectrum given an observation point, and energy range

  char const* SurfacePlane = "";
  size_t      NX1 = 0;
  size_t      NX2 = 0;
  double      Width_X1 = 0;
  double      Width_X2 = 0;
  PyObject*   List_NPoints= PyList_New(0);
  PyObject*   List_Width= PyList_New(0);
  PyObject*   List_Translation = PyList_New(0);
  PyObject*   List_Rotations = PyList_New(0);
  PyObject*   List_X0X1X2 = PyList_New(0);
  int         NormalDirection = 0;
  int         Dim = 2;
  double      Energy_eV = 0;
  char const* OutFileName = "";


  static char *kwlist[] = {"energy_eV", "npoints", "plane", "normal", "dim", "width", "rotations", "translation", "x0x1x2", "ofile", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "dO|siiOOOOs", kwlist,
                                                                &Energy_eV,
                                                                &List_NPoints,
                                                                &SurfacePlane,
                                                                &NormalDirection,
                                                                &Dim,
                                                                &List_Width,
                                                                &List_Rotations,
                                                                &List_Translation,
                                                                &List_X0X1X2,
                                                                &OutFileName)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Check requested dimension
  if (Dim != 2 & Dim != 3) {
    PyErr_SetString(PyExc_ValueError, "'dim' must be 2 or 3");
    return NULL;
  }


  // The rectangular surface object we'll use
  TSurfacePoints_Rectangle Surface;

  if (PyList_Size(List_NPoints) == 2) {
    // NPoints in [m]
    NX1 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 0));
    NX2 = PyInt_AsSsize_t(PyList_GetItem(List_NPoints, 1));
  } else {
    PyErr_SetString(PyExc_ValueError, "'npoints' must be [int, int]");
    return NULL;
  }



  if (NX1 <= 0 || NX2 <= 0) {
    PyErr_SetString(PyExc_ValueError, "an entry in 'npoints' is <= 0");
    return NULL;
  }

  // Vectors for rotations and translations.  Default to 0
  TVector3D Rotations(0, 0, 0);
  TVector3D Translation(0, 0, 0);

  // Check for Rotations in the input
  if (PyList_Size(List_Rotations) != 0) {
    try {
      Rotations = SRS_ListAsTVector3D(List_Rotations);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'rotations'");
      return NULL;
    }
  }


  // Check for Translation in the input
  if (PyList_Size(List_Translation) != 0) {
    try {
      Translation = SRS_ListAsTVector3D(List_Translation);
    } catch (std::length_error e) {
      PyErr_SetString(PyExc_ValueError, "Incorrect format in 'translation'");
      return NULL;
    }
  }

  if (PyList_Size(List_Width) == 2) {
    // Width in [m]
    Width_X1 = PyFloat_AsDouble(PyList_GetItem(List_Width, 0));
    Width_X2 = PyFloat_AsDouble(PyList_GetItem(List_Width, 1));
  }



  // If you are requesting a simple surface plane, check that you have widths
  if (SurfacePlane != "" && Width_X1 > 0 && Width_X2 > 0) {
    try {
      Surface.Init(SurfacePlane, NX1, NX2, Width_X1, Width_X2, Rotations, Translation, NormalDirection);
    } catch (std::invalid_argument e) {
      PyErr_SetString(PyExc_ValueError, e.what());
      return NULL;
    }
  }



  // If X0X1X2 defined
  std::vector<TVector3D> X0X1X2;

  if (PyList_Size(List_X0X1X2) != 0) {
    if (PyList_Size(List_X0X1X2) == 3) {
      for (int i = 0; i != 3; ++i) {
        PyObject* List_X = PyList_GetItem(List_X0X1X2, i);

        try {
          X0X1X2.push_back(SRS_ListAsTVector3D(List_X));
        } catch (std::length_error e) {
          PyErr_SetString(PyExc_ValueError, "Incorrect format in 'x0x1x2'");
          return NULL;
        }
      }
    } else {
      PyErr_SetString(PyExc_ValueError, "'x0x1x2' must have 3 XYZ points defined correctly");
      return NULL;
    }

    for (std::vector<TVector3D>::iterator it = X0X1X2.begin(); it != X0X1X2.end(); ++it) {
      it->RotateSelfXYZ(Rotations);
      *it += Translation;
    }

    // UPDATE: Check for orthogonality
    Surface.Init(NX1, NX2, X0X1X2[0], X0X1X2[1], X0X1X2[2], NormalDirection);
  }



  // Container for Point plus scalar
  T3DScalarContainer FluxContainer;


  // Actually calculate the spectrum
  bool const Directional = NormalDirection == 0 ? false : true;
  //self->obj->CalculateFlux(Surface, Energy_eV, FluxContainer, Dim, Directional);
  self->obj->CalculateFlux(Surface, Energy_eV, FluxContainer);


  // Build the output list of: [[[x, y, z], Flux], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = FluxContainer.GetNPoints();

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = FluxContainer.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    PyList_Append(PList2, SRS_TVector3DAsList(P.GetX()));
    PyList_Append(PList2, Py_BuildValue("f", P.GetV()));
    PyList_Append(PList, PList2);

  }

  return PList;
}




static PyObject* SRS_CalculateElectricFieldTimeDomain (SRSObject* self, PyObject* args, PyObject *keywds)
{
  // Calculate the electric field in the proper time domain.
  // The warning for using this function is that it returns unevenly spaced
  // time steps in the lab frame.  The calculation is based on the time stamps/steps
  // in the Trajectory object.  ie don't blindly employ a DFT or FFT.

  PyObject*   List_Obs = PyList_New(0);
  char const* OutFileName = "";


  static char *kwlist[] = {"obs", "ofile", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|O", kwlist,
                                                        &List_Obs,
                                                        &OutFileName)) {
    return NULL;
  }

  // Check if a beam is at least defined
  if (self->obj->GetNParticleBeams() < 1) {
    PyErr_SetString(PyExc_ValueError, "No particle beam defined");
    return NULL;
  }


  // Observation point
  TVector3D Obs(0, 0, 0);
  try {
      Obs = SRS_ListAsTVector3D(List_Obs);
  } catch (std::length_error e) {
    PyErr_SetString(PyExc_ValueError, "Incorrect format in 'obs'");
    return NULL;
  }

  T3DScalarContainer XYZT;
  self->obj->CalculateElectricFieldTimeDomain(Obs, XYZT);


  // Build the output list of: [[[x, y, z], Flux], [...]]
  // Create a python list
  PyObject *PList = PyList_New(0);

  size_t const NPoints = XYZT.GetNPoints();

  for (size_t i = 0; i != NPoints; ++i) {
    T3DScalar P = XYZT.GetPoint(i);

    // Inner list for each point
    PyObject *PList2 = PyList_New(0);


    // Add position and value to list
    PyList_Append(PList2, Py_BuildValue("f", P.GetV()));
    PyList_Append(PList2, SRS_TVector3DAsList(P.GetX()));
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

  {"pi",                                (PyCFunction) SRS_Pi,                              METH_NOARGS,  "do you want some pi?"},
  {"qe",                                (PyCFunction) SRS_Qe,                              METH_NOARGS,  "elementary charge in [C]"},
  {"me",                                (PyCFunction) SRS_Me,                              METH_NOARGS,  "electron mass in [kg]"},


  {"set_ctstart",                       (PyCFunction) SRS_SetCTStart,                      METH_O,       "set the start time in [m]"},
  {"get_ctstart",                       (PyCFunction) SRS_GetCTStart,                      METH_NOARGS,  "get the start time in [m]"},
  {"set_ctstop",                        (PyCFunction) SRS_SetCTStop,                       METH_O,       "set the stop time in [m]"},
  {"get_ctstop",                        (PyCFunction) SRS_GetCTStop,                       METH_NOARGS,  "get the stop time in [m]"},
  {"set_ctstartstop",                   (PyCFunction) SRS_SetCTStartStop,                  METH_VARARGS, "set the start and stop time in [m]"},
  {"set_npoints_trajectory",            (PyCFunction) SRS_SetNPointsTrajectory,            METH_O,       "set the total number of points for the trajectory"},
  {"get_npoints_trajectory",            (PyCFunction) SRS_GetNPointsTrajectory,            METH_NOARGS,  "get the total number of points for the trajectory"},
                                                                                          
  {"add_magnetic_field",                (PyCFunction) SRS_AddMagneticField,                METH_VARARGS | METH_KEYWORDS, "add a magnetic field from a file"},
  {"add_magnetic_field_function",       (PyCFunction) SRS_AddMagneticFieldFunction,        METH_VARARGS, "add a magnetic field in form of python function"},
  {"add_magnetic_field_gaussian",       (PyCFunction) SRS_AddMagneticFieldGaussian,        METH_VARARGS | METH_KEYWORDS, "add a magnetic field in form of 3D gaussian"},
  {"add_magnetic_field_uniform",        (PyCFunction) SRS_AddMagneticFieldUniform,         METH_VARARGS | METH_KEYWORDS, "add a uniform magnetic field in 3D"},
  {"add_undulator",                     (PyCFunction) SRS_AddMagneticFieldIdealUndulator,  METH_VARARGS | METH_KEYWORDS, "add magnetic field from ideal undulator in 3D"},
  {"get_bfield",                        (PyCFunction) SRS_GetBField,                       METH_VARARGS, "get the magnetic field at a given position in space (and someday time?)"},
  {"clear_magnetic_fields",             (PyCFunction) SRS_ClearMagneticFields,             METH_NOARGS,  "clear all internal magnetic fields"},
                                                                                          
                                                                                          
  {"set_particle_beam",                 (PyCFunction) SRS_SetParticleBeam,                 METH_VARARGS | METH_KEYWORDS, "add a particle beam"},
  {"add_particle_beam",                 (PyCFunction) SRS_AddParticleBeam,                 METH_VARARGS | METH_KEYWORDS, "add a particle beam"},
  {"clear_particle_beams",              (PyCFunction) SRS_ClearParticleBeams,              METH_NOARGS,  "Clear all existing particle beams from SRS"},
                                                                                          
  {"set_new_particle",                  (PyCFunction) SRS_SetNewParticle,                  METH_VARARGS | METH_KEYWORDS,  "Set the internal particle to a new random particle"},
  {"get_particle_x0",                   (PyCFunction) SRS_GetParticleX0,                   METH_NOARGS,  "Get the position at t0"},
  {"get_particle_beta0",                (PyCFunction) SRS_GetParticleBeta0,                METH_NOARGS,  "Get the beta at t0"},
  {"get_particle_e0",                   (PyCFunction) SRS_GetParticleE0,                   METH_NOARGS,  "Get the Energy at t0"},
                                                                                          
  {"calculate_trajectory",              (PyCFunction) SRS_CalculateTrajectory,             METH_NOARGS,  "Calclate the trajectory for the current particle"},
  {"get_trajectory",                    (PyCFunction) SRS_GetTrajectory,                   METH_NOARGS,  "Get the trajectory for the current particle as 2 3D lists [[x, y, z], [BetaX, BetaY, BetaZ]]"},

  {"calculate_spectrum",                (PyCFunction) SRS_CalculateSpectrum,               METH_VARARGS | METH_KEYWORDS, "calculate the spectrum at an observation point"},
  {"calculate_spectrum_gpu",            (PyCFunction) SRS_CalculateSpectrumGPU,            METH_VARARGS | METH_KEYWORDS, "calculate the spectrum at an observation point"},

  {"calculate_total_power",             (PyCFunction) SRS_CalculateTotalPower,             METH_NOARGS,  "calculate total power radiated"},
  {"calculate_power_density_rectangle", (PyCFunction) SRS_CalculatePowerDensityRectangle,  METH_VARARGS | METH_KEYWORDS, "calculate the power density given a surface"},
  {"calculate_power_density_rectangle_gpu", (PyCFunction) SRS_CalculatePowerDensityRectangleGPU,  METH_VARARGS | METH_KEYWORDS, "calculate the power density given a surface"},
  {"calculate_power_density",           (PyCFunction) SRS_CalculatePowerDensity,           METH_VARARGS | METH_KEYWORDS, "calculate the power density given a surface"},
  {"calculate_flux_rectangle",          (PyCFunction) SRS_CalculateFluxRectangle,          METH_VARARGS | METH_KEYWORDS, "calculate the flux given a surface"},

  {"calculate_electric_field",          (PyCFunction) SRS_CalculateElectricFieldTimeDomain,METH_VARARGS | METH_KEYWORDS, "calculate the electric field in the time domain"},

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


