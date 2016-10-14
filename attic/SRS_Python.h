#ifndef GUARD_SRS_Python_h
#define GUARD_SRS_Python_h
////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Jul 12 08:36:18 EDT 2016
//
// This is a header file for the SRS_Python 'SRS' module
//
////////////////////////////////////////////////////////////////////

// Include Python.h first!
#include <Python.h>

#include "SRS.h"

// The python SRS object
typedef struct {
  // Define the SRSObject struct which contains the class I want
  PyObject_HEAD
  SRS* obj;
} SRSObject;



static void SRS_dealloc(SRSObject* self);
static PyObject* SRS_new (PyTypeObject* type, PyObject* args, PyObject* kwds);
static TVector3D SRS_ListAsTVector3D (PyObject* List);
static PyObject* SRS_TVector3DAsList (TVector3D const& V);
static PyObject* SRS_Pi (SRSObject* self, PyObject* arg);
static PyObject* SRS_SetCTStart (SRSObject* self, PyObject* arg);
static PyObject* SRS_GetCTStart (SRSObject* self);
static PyObject* SRS_SetCTStop (SRSObject* self, PyObject* arg);
static PyObject* SRS_SetCTStop2 (SRSObject* self, PyObject* args, PyObject *keywds);
static PyObject* SRS_GetCTStop (SRSObject* self);
static PyObject* SRS_SetCTStartStop (SRSObject* self, PyObject* args);
static PyObject* SRS_GetNPointsTrajectory (SRSObject* self);
static PyObject* SRS_SetNPointsTrajectory (SRSObject* self, PyObject* arg);
static PyObject* SRS_AddMagneticField (SRSObject* self, PyObject* args, PyObject* keywds);
static PyObject* SRS_AddMagneticFieldFunction (SRSObject* self, PyObject* args);
static PyObject* SRS_AddMagneticFieldGaussian (SRSObject* self, PyObject* args, PyObject* keywds);
static PyObject* SRS_ClearMagneticFields (SRSObject* self);
static PyObject* SRS_GetBField (SRSObject* self, PyObject* args);
static PyObject* SRS_AddParticleBeam (SRSObject* self, PyObject* args, PyObject* keywds);
static PyObject* SRS_SetNewParticle (SRSObject* self, PyObject* args, PyObject* keywds);
static PyObject* SRS_ClearParticleBeams (SRSObject* self);
static PyObject* SRS_CalculateTrajectory (SRSObject* self);
static PyObject* SRS_GetTrajectory (SRSObject* self);
static PyObject* SRS_GetSpectrum (SRSObject* self);
static PyObject* SRS_CalculateSpectrum (SRSObject* self, PyObject* args, PyObject* keywds);
static PyObject* SRS_CalculateTotalPower (SRSObject* self);
static PyObject* SRS_CalculatePowerDensityRectangle (SRSObject* self, PyObject* args, PyObject *keywds);
static PyObject* SRS_CalculateFluxRectangle (SRSObject* self, PyObject* args, PyObject *keywds);
static PyObject* SRS_CalculateFluxRectangle2 (SRSObject* self, PyObject* args);























#endif
