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

#include "pylith/faults/KinSrcConstRate.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcConstRate::KinSrcConstRate(void) {
    pylith::utils::PyreComponent::setName("kinsrcconstrate");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcConstRate::~KinSrcConstRate(void) {}


// ------------------------------------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcConstRate::slipFn(const pylith::integer dim,
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
    const pylith::integer i_slipRate = 1;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* slipRate = &a[aOff[i_slipRate]];

    const pylith::integer i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        for (pylith::integer i = 0; i < dim; ++i) {
            slip[i] = slipRate[i] * (t - t0);
        } // for
    } else {
        for (pylith::integer i = 0; i < dim; ++i) {
            slip[i] = 0.0;
        } // for
    } // if/else

} // slipFn


// ------------------------------------------------------------------------------------------------
// Slip rate time function kernel.
void
pylith::faults::KinSrcConstRate::slipRateFn(const pylith::integer dim,
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
                                            PylithScalar slipRate[]) {
    const pylith::integer _numA = 2;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipRate);

    const pylith::integer i_initiationTime = 0;
    const pylith::integer i_slipRate = 1;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* slipRateAux = &a[aOff[i_slipRate]];

    const pylith::integer i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        for (pylith::integer i = 0; i < dim; ++i) {
            slipRate[i] = slipRateAux[i];
        } // for
    } else {
        for (pylith::integer i = 0; i < dim; ++i) {
            slipRate[i] = 0.0;
        } // for
    } // if/else

} // slipRateFn


// ------------------------------------------------------------------------------------------------
// Slip acceleration time function kernel.
void
pylith::faults::KinSrcConstRate::slipAccFn(const pylith::integer dim,
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
                                           PylithScalar slipAcc[]) {
    const pylith::integer _numA = 2;

    assert(_numA == numA);
    assert(slipAcc);

    for (pylith::integer i = 0; i < dim; ++i) {
        slipAcc[i] = 0.0;
    } // for

} // slipRateFn


// ------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcConstRate::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                      const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time function
    // kernel.

    _auxiliaryFactory->addInitiationTime(); // 0
    _auxiliaryFactory->addSlipRate(); // 1

    _slipFnKernel = pylith::faults::KinSrcConstRate::slipFn;
    _slipRateFnKernel = pylith::faults::KinSrcConstRate::slipRateFn;
    _slipAccFnKernel = pylith::faults::KinSrcConstRate::slipAccFn;

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// End of file
