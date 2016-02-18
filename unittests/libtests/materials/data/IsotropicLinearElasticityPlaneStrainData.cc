// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include "IsotropicLinearElasticityPlaneStrainData.hh"


const int pylith::materials::IsotropicLinearElasticityPlaneStrainData::numKernelsResidual = 2;
const int pylith::materials::IsotropicLinearElasticityPlaneStrainData::numKernelsJacobian = 4;

// ----------------------------------------------------------------------
// Constructor
pylith::materials::IsotropicLinearElasticityPlaneStrainData::IsotropicLinearElasticityPlaneStrainData(void) :
  filenameMesh(0),
  label(NULL),
  id(0),
  dimension(0),

  useInertia(false),
  useBodyForce(false),

  numSolnFields(0),
  discretizations(NULL),

  filenameAuxFieldsDB(NULL),

  kernelsRHSResidual(NULL),
  kernelsRHSJacobian(NULL),
  kernelsLHSResidual(NULL),
  kernelsLHSJacobianImplicit(NULL),
  kernelsLHSJacobianExplicit(NULL),

  querySolutionDisplacement(NULL),
  querySolutionVelocity(NULL),

  lengthScale(0),
  timeScale(0),
  pressureScale(0),
  densityScale(0)
{ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::materials::IsotropicLinearElasticityPlaneStrainData::~IsotropicLinearElasticityPlaneStrainData(void)
{ // destructor
} // destructor

// End of file