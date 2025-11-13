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

#include "pylith/feassemble/JacobianValues.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/DSLabelAccess.hh" // USES DSLabelAccess

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

namespace pylith {
    namespace feassemble {
        class _JacobianValues {
public:

            static
            void evaluateKernel(scalar_array* cellMat,
                                const JacobianValues::JacobianKernel& kernel,
                                const pylith::real t,
                                const pylith::real dt,
                                const pylith::real s_tshift,
                                const pylith::topology::Field& solution,
                                const pylith::feassemble::DSLabelAccess& dsLabel);

        };
    }
}

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::JacobianValues::JacobianValues(void) {
    GenericComponent::setName("jacobianvalues");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::JacobianValues::~JacobianValues(void) {}


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::JacobianValues::setKernels(const std::vector<JacobianKernel>& kernelsJacobian,
                                               const std::vector<JacobianKernel>& kernelsPrecond) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setKernels(# Jacobian kernels="<<kernelsJacobian.size()<<", # preconditioner kernels="<<kernelsPrecond.size()<<")");

    _kernelsJacobian = kernelsJacobian;
    _kernelsPrecond = kernelsPrecond;

    PYLITH_METHOD_END;
} // setKernels


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::JacobianValues::computeLHSJacobian(PetscMat jacobianMat,
                                                       PetscMat precondMat,
                                                       const pylith::real t,
                                                       const pylith::real dt,
                                                       const pylith::real s_tshift,
                                                       const pylith::topology::Field& solution,
                                                       const pylith::feassemble::DSLabelAccess& dsLabel) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("computeLHSJacobian (jacobianMat = "<<jacobianMat<<", precondMat = "<<precondMat<<", t = "<<t<<", dt = "<<dt<<", solution = "<<solution.getName()<<") ");

    PetscDM dm = dsLabel.dm();assert(dm);

    PetscErrorCode err = PETSC_SUCCESS;
    const pylith::integer numCells = dsLabel.numCells();
    const pylith::integer* cellIndices = nullptr;
    err = ISGetIndices(dsLabel.cellsIS(), &cellIndices);PYLITH_CHECK_ERROR(err);

    pylith::integer totalDof = 0;
    err = PetscDSGetTotalDimension(dsLabel.ds(), &totalDof);PYLITH_CHECK_ERROR(err);
    pylith::scalar_array cellMat(totalDof*totalDof);

    if (jacobianMat) {
        for (pylith::integer iCell = 0; iCell < numCells; ++iCell) {
            const pylith::integer cell = cellIndices[iCell];
            cellMat = 0.0;

            for (size_t i = 0; i < _kernelsJacobian.size(); ++i) {
                _JacobianValues::evaluateKernel(&cellMat, _kernelsJacobian[i], t, dt, s_tshift, solution, dsLabel);
            } // for
            err = DMPlexMatSetClosure(dsLabel.dm(), nullptr, nullptr, jacobianMat, cell, &cellMat[0],
                                      INSERT_VALUES);PYLITH_CHECK_ERROR(err);
        } // for
    } // if

    if (precondMat && (jacobianMat != precondMat)) {
        for (pylith::integer iCell = 0; iCell < numCells; ++iCell) {
            const pylith::integer cell = cellIndices[iCell];
            cellMat = 0.0;

            for (size_t i = 0; i < _kernelsPrecond.size(); ++i) {
                _JacobianValues::evaluateKernel(&cellMat, _kernelsPrecond[i], t, dt, s_tshift, solution, dsLabel);
            } // for
            err = DMPlexMatSetClosure(dsLabel.dm(), nullptr, nullptr, precondMat, cell, &cellMat[0],
                                      INSERT_VALUES);PYLITH_CHECK_ERROR(err);
        } // for
    } // if
    err = ISRestoreIndices(dsLabel.cellsIS(), &cellIndices);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // computeLHSJacobian


// ------------------------------------------------------------------------------------------------
// Compute Jacobian values associated with s_tshift on diagonal.
void
pylith::feassemble::JacobianValues::blockDiag_tshift(pylith::scalar_array* cellMat,
                                                     const pylith::real t,
                                                     const pylith::real dt,
                                                     const pylith::real s_tshift,
                                                     const pylith::integer trialDof,
                                                     const pylith::integer trialOff,
                                                     const pylith::integer basisDof,
                                                     const pylith::integer basisOff,
                                                     const pylith::integer totalDim) {
    assert(trialDof == basisDof);

    const pylith::real value = s_tshift;
    for (pylith::integer iDof = 0; iDof < trialDof; ++iDof) {
        const pylith::integer index = (trialOff+iDof) * totalDim + basisOff + iDof;
        (*cellMat)[index] = value;
    } // for
} // blockDiag_tshift


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::_JacobianValues::evaluateKernel(scalar_array* cellMat,
                                                    const JacobianValues::JacobianKernel& kernel,
                                                    const pylith::real t,
                                                    const pylith::real dt,
                                                    const pylith::real s_tshift,
                                                    const pylith::topology::Field& solution,
                                                    const pylith::feassemble::DSLabelAccess& dsLabel) {
    const pylith::topology::Field::SubfieldInfo& infoTrial = solution.getSubfieldInfo(kernel.subfieldTrial.c_str());
    const pylith::topology::Field::SubfieldInfo& infoBasis = solution.getSubfieldInfo(kernel.subfieldBasis.c_str());
    const size_t i_trial = infoTrial.index;
    const size_t i_basis = infoBasis.index;

    PetscErrorCode err = PETSC_SUCCESS;
    pylith::integer trialOff, trialDof, basisOff, basisDof;
    PetscFE fe = nullptr;
    err = PetscDSGetFieldOffset(dsLabel.ds(), i_trial, &trialOff);PYLITH_CHECK_ERROR(err);
    err = PetscDSGetDiscretization(dsLabel.ds(), i_trial, (PetscObject*) &fe);PYLITH_CHECK_ERROR(err);
    err = PetscFEGetDimension(fe, &trialDof);PYLITH_CHECK_ERROR(err);

    err = PetscDSGetFieldOffset(dsLabel.ds(), i_basis, &basisOff);PYLITH_CHECK_ERROR(err);
    err = PetscDSGetDiscretization(dsLabel.ds(), i_basis, (PetscObject*) &fe);PYLITH_CHECK_ERROR(err);
    err = PetscFEGetDimension(fe, &basisDof);PYLITH_CHECK_ERROR(err);

    pylith::integer totalDof = 0;
    err = PetscDSGetTotalDimension(dsLabel.ds(), &totalDof);PYLITH_CHECK_ERROR(err);

    kernel.function(cellMat, t, dt, s_tshift, trialDof, trialOff, basisDof, basisOff, totalDof);
} // evaluateKernel


// End of file
