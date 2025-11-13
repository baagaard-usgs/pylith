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

#include "pylith/feassemble/ConstraintCxxFn.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()
#include <stdexcept> // USES std::runtime_error

namespace pylith {
    namespace feassemble {
        // Trampoline class for ConstraintCxxFn::create() factory method
        class ConstraintCxxFnWrap : public ConstraintCxxFn {
public:

            ConstraintCxxFnWrap(const std::shared_ptr<pylith::problems::Physics>& physics) : ConstraintCxxFn(physics) {}


        };

        class _ConstraintCxxFn {
public:

            static
            void setSolution(const pylith::topology::Field* field,
                             const pylith::real t,
                             PetscUserFieldFunc fn,
                             const pylith::feassemble::ConstraintCxxFn& constraint);

        };
    }
}

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::ConstraintCxxFn::ConstraintCxxFn(const std::shared_ptr<pylith::problems::Physics>& physics) :
    Constraint(physics),
    _fn(nullptr),
    _fnDot(nullptr) {
    GenericComponent::setName("constraintuserfn");
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::feassemble::ConstraintCxxFn>
pylith::feassemble::ConstraintCxxFn::create(const std::shared_ptr<pylith::problems::Physics>& physics) {
    return std::make_shared<ConstraintCxxFnWrap>(physics);
} // create


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::ConstraintCxxFn::~ConstraintCxxFn(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Set constraint kernel.
void
pylith::feassemble::ConstraintCxxFn::setCxxFn(const PetscUserFieldFunc fn) {
    _fn = fn;
} // setCxxFn


// ------------------------------------------------------------------------------------------------
// Set constraint kernel time derivative.
void
pylith::feassemble::ConstraintCxxFn::setCxxFnDot(const PetscUserFieldFunc fnDot) {
    _fnDot = fnDot;
} // setCxxFnDot


// ------------------------------------------------------------------------------------------------
// Initialize constraint domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::ConstraintCxxFn::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" initialize(solution="<<solution.getName()<<")");

    Constraint::initialize(solution);

    const pylith::problems::Observer::NotificationType notification = pylith::problems::ObserverPhysics::DIAGNOSTIC;
    _observers->notifyObservers(0.0, 0, solution, notification);

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a nullptr
    // label) is the correct one.
    PetscErrorCode err = PETSC_SUCCESS;
    PetscDS prob = nullptr;
    DMLabel label = nullptr;
    void* context = nullptr;
    err = DMGetDS(solution.getDM(), &prob);PYLITH_CHECK_ERROR(err);
    const pylith::integer i_field = solution.getSubfieldInfo(_subfieldName.c_str()).index;
    err = DMGetLabel(solution.getDM(), _labelName.c_str(), &label);PYLITH_CHECK_ERROR(err);
    err = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, _labelName.c_str(), label, 1, &_labelValue, i_field,
                             _constrainedDOF.size(), &_constrainedDOF[0], (void (*)(void)) _fn, (void (*)(void)) _fnDot, context, nullptr);
    PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::ConstraintCxxFn::setSolution(pylith::feassemble::IntegrationData& integrationData) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setSolution(integrationData="<<integrationData.str()<<")");

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::real t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);

    _ConstraintCxxFn::setSolution(solution, t, _fn, *this);

    if (_fnDot && integrationData.hasField(pylith::feassemble::IntegrationData::solution_dot)) {
        const pylith::topology::Field* solutionDot = integrationData.getField(pylith::feassemble::IntegrationData::solution_dot);
        assert(solutionDot);
        _ConstraintCxxFn::setSolution(solutionDot, t, _fnDot, *this);
    } // if

    PYLITH_METHOD_END;
} // setSolution


// ------------------------------------------------------------------------------------------------
// Set constrained values in solution field.
void
pylith::feassemble::_ConstraintCxxFn::setSolution(const pylith::topology::Field* field,
                                                  const pylith::real t,
                                                  PetscUserFieldFunc fn,
                                                  const pylith::feassemble::ConstraintCxxFn& constraint) {
    PYLITH_METHOD_BEGIN;
    assert(field);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmField = field->getDM();

    // Get label for constraint.
    PetscDMLabel dmLabel = nullptr;
    err = DMGetLabel(dmField, constraint._labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);

    void* context = nullptr;
    const int fieldIndex = field->getSubfieldInfo(constraint._subfieldName.c_str()).index;
    const pylith::integer numConstrained = constraint._constrainedDOF.size();
    assert(field->getLocalVector());
    err = DMPlexLabelAddCells(dmField, dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexInsertBoundaryValuesEssential(dmField, t, fieldIndex, numConstrained, &constraint._constrainedDOF[0], dmLabel, 1,
                                              &constraint._labelValue, fn, context, field->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelClearCells(dmField, dmLabel);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // setSolution


// End of file
