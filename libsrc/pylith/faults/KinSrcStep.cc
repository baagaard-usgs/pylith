// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/faults/KinSrcStep.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcStep::KinSrcStep(void) {
    pylith::utils::PyreComponent::setName("kinsrcstep");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcStep::~KinSrcStep(void) {}


#include <iostream>
// ---------------------------------------------------------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcStep::slipFn(const pylith::integer dim,
                                   const pylith::integer numS,
                                   const pylith::integer numA,
                                   const pylith::integer sOff[],
                                   const pylith::integer sOff_x[],
                                   const PylithScalar s[],
                                   const PylithScalar s_t[],
                                   const PylithScalar s_x[],
                                   const pylith::integer aOff[],
                                   const pylith::integer aOff_x[],
                                   const PylithScalar a[],
                                   const PylithScalar a_t[],
                                   const PylithScalar a_x[],
                                   const pylith::real t,
                                   const PylithScalar x[],
                                   const pylith::integer numConstants,
                                   const PylithScalar constants[],
                                   PylithScalar slip[]) {
    const pylith::integer _numA = 2;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slip);

    const pylith::integer i_initiationTime = 0;
    const pylith::integer i_finalSlip = 1;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];

    const pylith::integer i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        for (pylith::integer i = 0; i < dim; ++i) {
            slip[i] = finalSlip[i];
        } // for
    } else {
        for (pylith::integer i = 0; i < dim; ++i) {
            slip[i] = 0.0;
        } // for
    } // for

} // slipFn


// ---------------------------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcStep::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                 const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time function
    // kernel.

    _auxiliaryFactory->addInitiationTime(); // 0
    _auxiliaryFactory->addFinalSlip(); // 1

    _slipFnKernel = pylith::faults::KinSrcStep::slipFn;
    _slipRateFnKernel = nullptr; // Undefined for step function.
    _slipAccFnKernel = nullptr; // Undefined for step function.

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// End of file
