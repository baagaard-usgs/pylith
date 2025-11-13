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

#include "pylith/feassemble/ConstraintSimple.hh" // implementation of object methods

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

namespace pylith {
    namespace feassemble {
        // Trampoline class for ConstraintSimple::create() factory method
        class ConstraintSimpleWrap : public ConstraintSimple {
public:

            ConstraintSimpleWrap(const std::shared_ptr<pylith::problems::Physics>& physics) : ConstraintSimple(physics) {}


        };
    } // feassemble
} // pylith

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintSimple::ConstraintSimple(const std::shared_ptr<pylith::problems::Physics>& physics) :
    Constraint(physics),
    _fn(nullptr) {
    GenericComponent::setName("constraintsimple");
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::feassemble::ConstraintSimple>
pylith::feassemble::ConstraintSimple::create(const std::shared_ptr<pylith::problems::Physics>& physics) {
    return std::make_shared<ConstraintSimpleWrap>(physics);
} // create


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintSimple::~ConstraintSimple(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Set constraint kernel.
void
pylith::feassemble::ConstraintSimple::setUserFn(const PetscUserFieldFunc fn) {
    _fn = fn;
} // setSimple


// ------------------------------------------------------------------------------------------------
// Initialize constraint domain. Update observers.
void
pylith::feassemble::ConstraintSimple::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(solution="<<solution.getName()<<")");

    assert(_physics);
    _observers = nullptr;

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dm = solution.getDM();
    PetscDMLabel label = nullptr;
    PetscDS ds = nullptr;
    void* context = nullptr;
    pylith::integer i_field = -1;
    pylith::integer *closure = nullptr;
    PetscIS pointIS;
    const pylith::integer* points = nullptr;
    pylith::integer point, cStart, cEnd, clSize;
    err = DMGetLabel(dm, _labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumIS(label, _labelValue, &pointIS);PYLITH_CHECK_ERROR(err);
    if (!pointIS) {
        PYLITH_METHOD_END;
    } // if
    err = ISGetIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);assert(points);
    point = points[0];
    err = ISRestoreIndices(pointIS, &points);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&pointIS);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetTransitiveClosure(dm, point, PETSC_FALSE, &clSize, &closure);PYLITH_CHECK_ERROR(err);
    for (pylith::integer cl = 0; cl < clSize*2; cl += 2) {
        PetscDS cds;
        const pylith::integer q = closure[cl];
        pylith::integer Nf;

        if ((q < cStart) || (q >= cEnd)) { continue;}
        err = DMGetCellDS(dm, q, &cds, nullptr);PYLITH_CHECK_ERROR(err);
        err = PetscDSGetNumFields(cds, &Nf);PYLITH_CHECK_ERROR(err);
        for (int f = 0; f < Nf; ++f) {
            PetscObject disc;
            const char *name;

            err = PetscDSGetDiscretization(cds, f, &disc);PYLITH_CHECK_ERROR(err);
            err = PetscObjectGetName(disc, &name);PYLITH_CHECK_ERROR(err);
            if (_subfieldName == std::string(name)) {ds = cds;i_field = f;break;}
        } // for
    } // for
    pylith::integer numConstrainedDOF = _constrainedDOF.size();
    pylith::integer* constrainedDOF = &_constrainedDOF[0];
    if (!ds) {
        // :KLUDGE: It is possible for a process to have a DOF that we need to constrain, but the process
        // may not have any cells with that DOF. The underlying code doesn't actually care if the point is
        // in the DS, so just get any DS and use it for the constraint.
        err = DMGetDS(dm, &ds);PYLITH_CHECK_ERROR(err);
        pylith::integer numDS = 0, numFields = 0;
        i_field = solution.getSubfieldInfo(_subfieldName.c_str()).index;
        numConstrainedDOF = 0;
        constrainedDOF = nullptr;

        err = DMGetNumDS(dm, &numDS);PYLITH_CHECK_ERROR(err);
        for (pylith::integer s = 0; s < numDS; ++s) {
            err = DMGetRegionNumDS(dm, s, nullptr, nullptr, &ds, nullptr);PYLITH_CHECK_ERROR(err);
            err = PetscDSGetNumFields(ds, &numFields);PYLITH_CHECK_ERROR(err);
            if (i_field < numFields) { break;}
        } // for

    } // if
    err = DMPlexRestoreTransitiveClosure(solution.getDM(), point, PETSC_FALSE, &clSize, &closure);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(solution.getDM(), _labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscDSAddBoundary(ds, DM_BC_ESSENTIAL, _labelName.c_str(), label, 1, &_labelValue, i_field,
                             numConstrainedDOF, constrainedDOF, (void (*)(void)) _fn, nullptr, context, nullptr);
    PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dm, nullptr, "-constraint_simple_dm_view");PYLITH_CHECK_ERROR(err);
    {
        pylith::integer numDS;
        err = DMGetNumDS(dm, &numDS);PYLITH_CHECK_ERROR(err);
        for (int s = 0; s < numDS; ++s) {
            err = DMGetRegionNumDS(dm, s, nullptr, nullptr, &ds, nullptr);PYLITH_CHECK_ERROR(err);
            err = PetscObjectViewFromOptions((PetscObject) ds, nullptr, "-constraint_simple_ds_view");PYLITH_CHECK_ERROR(err);
        }
    }

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::ConstraintSimple::setSolution(pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSolution(integrationData="<<integrationData.str()<<")");

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::real t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmSoln = solution->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = nullptr;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    void* context = nullptr;
    const pylith::integer _labelValue = 1;
    const pylith::integer fieldIndex = solution->getSubfieldInfo(_subfieldName.c_str()).index;
    const pylith::integer numConstrained = _constrainedDOF.size();
    assert(solution->getLocalVector());
    err = DMPlexLabelAddCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssential(dmSoln, t, fieldIndex, numConstrained, &_constrainedDOF[0], dmLabel, 1,
                                              &_labelValue, _fn, context, solution->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmSoln, dmLabel);PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Displaying solution field");
        solution->view("solution field");
    } // if

    PYLITH_METHOD_END;
} // setSolution


// End of file
