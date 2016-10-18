#ifndef GUARD_OSCARS_Python_h
#define GUARD_OSCARS_Python_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 12 08:36:18 EDT 2016
//
// This is a header file for the OSCARS_Python 'OSCARS' module
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "OSCARS.h"

// The python OSCARS object
typedef struct {
  // Define the OSCARSObject struct which contains the class I want
  PyObject_HEAD
  OSCARS* obj;
} OSCARSObject;



static void OSCARS_dealloc(OSCARSObject* self);
static PyObject* OSCARS_new (PyTypeObject* type, PyObject* args, PyObject* kwds);
static TVector3D OSCARS_ListAsTVector3D (PyObject* List);
static PyObject* OSCARS_TVector3DAsList (TVector3D const& V);
static PyObject* OSCARS_Pi (OSCARSObject* self, PyObject* arg);
static PyObject* OSCARS_GetCTStart (OSCARSObject* self);
static PyObject* OSCARS_GetCTStop (OSCARSObject* self);
static PyObject* OSCARS_SetCTStartStop (OSCARSObject* self, PyObject* args);
static PyObject* OSCARS_GetNPointsTrajectory (OSCARSObject* self);
static PyObject* OSCARS_SetNPointsTrajectory (OSCARSObject* self, PyObject* arg);
static PyObject* OSCARS_AddMagneticField (OSCARSObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARS_AddMagneticFieldFunction (OSCARSObject* self, PyObject* args);
static PyObject* OSCARS_AddMagneticFieldGaussian (OSCARSObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARS_ClearMagneticFields (OSCARSObject* self);
static PyObject* OSCARS_GetBField (OSCARSObject* self, PyObject* args);
static PyObject* OSCARS_AddParticleBeam (OSCARSObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARS_SetNewParticle (OSCARSObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARS_ClearParticleBeams (OSCARSObject* self);
static PyObject* OSCARS_CalculateTrajectory (OSCARSObject* self);
static PyObject* OSCARS_GetTrajectory (OSCARSObject* self);
static PyObject* OSCARS_GetSpectrum (OSCARSObject* self);
static PyObject* OSCARS_CalculateSpectrum (OSCARSObject* self, PyObject* args, PyObject* keywds);
static PyObject* OSCARS_CalculateTotalPower (OSCARSObject* self);
static PyObject* OSCARS_CalculatePowerDensityRectangle (OSCARSObject* self, PyObject* args, PyObject *keywds);
static PyObject* OSCARS_CalculateFluxRectangle (OSCARSObject* self, PyObject* args, PyObject *keywds);























#endif
