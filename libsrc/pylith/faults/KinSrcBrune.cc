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

#include "pylith/faults/KinSrcBrune.hh" // implementation of object methods

#include "pylith/faults/KinSrcAuxiliaryFactory.hh" // USES KinSrcAuxiliaryFactory

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()

namespace pylith {
    namespace faults {
        namespace _KinSrcBrune {
            inline
            pylith::real
            tau(const pylith::real riseTime) {
                return 0.21081916 * riseTime;
            } // tau


        } // _KinSrcBrune
    } // faults
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::faults::KinSrcBrune::KinSrcBrune(void) {
    pylith::utils::PyreComponent::setName("kinsrcbrune");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::faults::KinSrcBrune::~KinSrcBrune(void) {}


// ------------------------------------------------------------------------------------------------
// Slip time function kernel.
void
pylith::faults::KinSrcBrune::slipFn(const pylith::integer dim,
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
    const pylith::integer _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slip);

    const pylith::integer i_initiationTime = 0;
    const pylith::integer i_finalSlip = 1;
    const pylith::integer i_riseTime = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];

    const pylith::integer i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        const PylithScalar tau = _KinSrcBrune::tau(riseTime);
        for (pylith::integer i = 0; i < dim; ++i) {
            slip[i] = finalSlip[i] * (1.0 - exp(-(t-t0)/tau) * (1.0 + (t-t0)/tau));
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
pylith::faults::KinSrcBrune::slipRateFn(const pylith::integer dim,
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
    const pylith::integer _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipRate);

    const pylith::integer i_initiationTime = 0;
    const pylith::integer i_finalSlip = 1;
    const pylith::integer i_riseTime = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];

    const pylith::integer i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        const PylithScalar tau = _KinSrcBrune::tau(riseTime);
        for (pylith::integer i = 0; i < dim; ++i) {
            slipRate[i] = finalSlip[i] * (t-t0)/(tau*tau) * exp(-(t-t0)/tau);
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
pylith::faults::KinSrcBrune::slipAccFn(const pylith::integer dim,
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
    const pylith::integer _numA = 3;

    assert(_numA == numA);
    assert(aOff);
    assert(a);
    assert(slipAcc);

    const pylith::integer i_initiationTime = 0;
    const pylith::integer i_finalSlip = 1;
    const pylith::integer i_riseTime = 2;
    const PylithScalar initiationTime = a[aOff[i_initiationTime]];
    const PylithScalar* finalSlip = &a[aOff[i_finalSlip]];
    const PylithScalar riseTime = a[aOff[i_riseTime]];

    const pylith::integer i_originTime = 0;
    const PylithScalar originTime = constants[i_originTime];
    const PylithScalar t0 = originTime + initiationTime;

    if (t >= t0) {
        const PylithScalar tau = _KinSrcBrune::tau(riseTime);
        for (pylith::integer i = 0; i < dim; ++i) {
            slipAcc[i] = finalSlip[i] * 1.0/(tau*tau) * (1.0 - (t-t0)/tau) * exp(-(t-t0)/tau);
        } // for
    } else {
        for (pylith::integer i = 0; i < dim; ++i) {
            slipAcc[i] = 0.0;
        } // for
    } // if/else
} // slipAccFn


// ------------------------------------------------------------------------------------------------
// Preinitialize earthquake source. Set names/sizes of auxiliary subfields.
void
pylith::faults::KinSrcBrune::_auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                                                  const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_auxiliaryFieldSetup()");

    assert(_auxiliaryFactory);
    assert(cs);
    _auxiliaryFactory->initialize(_auxiliaryField, normalizer, cs->getSpaceDim());

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the slip time
    // function
    // kernel.

    _auxiliaryFactory->addInitiationTime(); // 0
    _auxiliaryFactory->addFinalSlip(); // 1
    _auxiliaryFactory->addRiseTime(); // 2

    _slipFnKernel = pylith::faults::KinSrcBrune::slipFn;
    _slipRateFnKernel = pylith::faults::KinSrcBrune::slipRateFn;

    PYLITH_METHOD_END;
} // _auxiliaryFieldSetup


// End of file
