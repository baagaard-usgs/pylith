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

#include "pylith/feassemble/Constraint.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/problems/ObserversPhysics.hh" // USES ObserversPhysics
#include "pylith/problems/Physics.hh" // USES Physics
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::Constraint::Constraint(const std::shared_ptr<pylith::problems::Physics>& physics) :
    PhysicsImplementation(physics),
    _subfieldName(""),
    _labelName(""),
    _labelValue(1) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::Constraint::~Constraint(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::Constraint::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _boundaryMesh.reset();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set name of constrained solution subfield.
void
pylith::feassemble::Constraint::setSubfieldName(const char* value) {
    PYLITH_JOURNAL_DEBUG("setSubfieldName(value="<<value<<")");

    if (!value || (0 == strlen(value))) {
        std::ostringstream msg;
        assert(_physics);
        msg << "Empty string given for name of solution subfield for constraint '" << _physics->getIdentifier()
            <<"'.";
        throw std::runtime_error(msg.str());
    } // if
    _subfieldName = value;
} // setSubfieldName


// ------------------------------------------------------------------------------------------------
// Get name of constrained solution subfield.
const char*
pylith::feassemble::Constraint::getSubfieldName(void) const {
    return _subfieldName.c_str();
} // getSubfieldName


// ------------------------------------------------------------------------------------------------
// Set name of label marking boundary associated with constraint.
void
pylith::feassemble::Constraint::setLabelName(const char* value) {
    PYLITH_JOURNAL_DEBUG("setLabelName(value="<<value<<")");

    if (strlen(value) == 0) {
        throw std::runtime_error("Empty string given for constraint label.");
    } // if

    _labelName = value;
} // setLabelName


// ------------------------------------------------------------------------------------------------
// Get name of label marking boundary associated with constraint.
const char*
pylith::feassemble::Constraint::getLabelName(void) const {
    return _labelName.c_str();
} // getLabelName


// ------------------------------------------------------------------------------------------------
// Set value of label marking boundary associated with constraint.
void
pylith::feassemble::Constraint::setLabelValue(const int value) {
    _labelValue = value;
} // setLabelValue


// ------------------------------------------------------------------------------------------------
// Get value of label marking boundary associated with constraint.
int
pylith::feassemble::Constraint::getLabelValue(void) const {
    return _labelValue;
} // getLabelValue


// ------------------------------------------------------------------------------------------------
// Set indices of constrained degrees of freedom at each location.
void
pylith::feassemble::Constraint::setConstrainedDOF(const pylith::integer_array& dof) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setConstrainedDOF(#dof="<<dof.size()<<")");

    const size_t size = dof.size();
    _constrainedDOF.resize(size);
    for (int i = 0; i < size; ++i) {
        if (dof[i] < 0) {
            std::ostringstream msg;
            assert(_physics);
            msg << "Constrained DOF '" << dof[i] << "' must be nonnegative in constraint component '" << _physics->getIdentifier() << "'.";
            throw std::runtime_error(msg.str());
        } // if
        _constrainedDOF[i] = dof[i];
    } // for

    PYLITH_METHOD_END;
} // setConstrainedDOF


// ------------------------------------------------------------------------------------------------
// Get indices of constrained degrees of freedom.
const pylith::integer_array&
pylith::feassemble::Constraint::getConstrainedDOF(void) const {
    return _constrainedDOF;
} // getConstrainedDOF


// ------------------------------------------------------------------------------------------------
// Get mesh associated with constrained boundary.
const pylith::topology::Mesh&
pylith::feassemble::Constraint::getPhysicsDomainMesh(void) const {
    assert(_boundaryMesh);
    return *_boundaryMesh;
} // getPhysicsDomainMesh


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::Constraint::setKernelsDiagnosticField(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsDiagnosticField(# kernels="<<kernels.size()<<")");

    _kernelsDiagnosticField = kernels;

    PYLITH_METHOD_END;
} // setKernelsDiagnosticField


// ----------------------------------------------------------------------
// Initialize boundary condition.
void
pylith::feassemble::Constraint::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(solution="<<solution.getName()<<")");

    const char* componentName = _physics->getFullIdentifier();
    _boundaryMesh = pylith::topology::MeshOps::createLowerDimMesh(solution.getMesh(), _labelName.c_str(), _labelValue, componentName);
    assert(_boundaryMesh);
    PetscDM dmBoundary = _boundaryMesh->getDM();assert(dmBoundary);
    pylith::topology::CoordsVisitor::optimizeClosure(dmBoundary);

    assert(_physics);
    _observers = _physics->getObservers();
    if (_observers) {
        _observers->setPhysicsImplementation(shared_from_this());
        _observers->setTimeScale(_physics->getNormalizer().getTimeScale());
    } // if

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Update at end of time step.
void
pylith::feassemble::Constraint::poststep(const pylith::real t,
                                         const pylith::integer tindex,
                                         const pylith::real dt,
                                         const pylith::topology::Field& solution,
                                         const pylith::problems::Observer::NotificationType notification) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("poststep(t="<<t<<", dt="<<dt<<")");

    notifyObservers(t, tindex, solution, notification);

    PYLITH_METHOD_END;
} // poststep


// ---------------------------------------------------------------------------------------------------------------------
// Set auxiliary field values for current time.
void
pylith::feassemble::Constraint::setState(const pylith::real t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setState(t="<<t<<") empty method");

    // Default is to do nothing.

    PYLITH_METHOD_END;
} // setState


// ------------------------------------------------------------------------------------------------
// Compute diagnostic field from auxiliary field.
void
pylith::feassemble::Constraint::_computeDiagnosticField(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeDiagnosticField()");

    if (!_diagnosticField || !_auxiliaryField) {
        PYLITH_METHOD_END;
    } // if

    assert(_auxiliaryField);
    assert(_diagnosticField);
    const pylith::real t = 0.0;
    const pylith::real dt = 0.0;
    _setKernelConstants(*_auxiliaryField, dt);

    const size_t numKernels = _kernelsDiagnosticField.size();
    assert(numKernels > 0);
    PetscBdPointFunc* kernelsArray = (numKernels > 0) ? new PetscBdPointFunc[numKernels] : nullptr;
    for (size_t iKernel = 0; iKernel < numKernels; ++iKernel) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _diagnosticField->getSubfieldInfo(_kernelsDiagnosticField[iKernel].subfield.c_str());
        kernelsArray[sinfo.index] = _kernelsDiagnosticField[iKernel].f;
    } // for

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM diagnosticDM = _diagnosticField->getDM();
    PetscDMLabel diagnosticFieldLabel = nullptr;
    const pylith::integer labelValue = 1;
    err = DMGetLabel(diagnosticDM, "output", &diagnosticFieldLabel);PYLITH_CHECK_ERROR(err);
    err = DMProjectBdFieldLabelLocal(diagnosticDM, t, diagnosticFieldLabel, 1, &labelValue, PETSC_DETERMINE, nullptr, _auxiliaryField->getLocalVector(), kernelsArray, INSERT_VALUES, _diagnosticField->getLocalVector());PYLITH_CHECK_ERROR(err);
    delete[] kernelsArray;kernelsArray = nullptr;

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing diagnostic field.");
        _diagnosticField->view("Diagnostic field");
    } // if

    PYLITH_METHOD_END;
} // _computeDiagnosticField


// ---------------------------------------------------------------------------------------------------------------------
// Set constants used in finite-element kernels (point-wise functions).
void
pylith::feassemble::Constraint::_setKernelConstants(const pylith::topology::Field& solution,
                                                    const pylith::real dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setKernelConstants(solution="<<solution.getName()<<", dt="<<dt<<")");

    assert(_physics);
    const pylith::real_array& constants = _physics->getKernelConstants(dt);

    PetscDS prob = nullptr;
    PetscDM dmSoln = solution.getDM();assert(dmSoln);

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a nullptr label)
    // is
    // the correct one.
    PetscErrorCode err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);assert(prob);
    if (constants.size() > 0) {
        err = PetscDSSetConstants(prob, constants.size(), const_cast<double*>(&constants[0]));PYLITH_CHECK_ERROR(err);
    } else {
        err = PetscDSSetConstants(prob, 0, nullptr);PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setKernelConstants


// End of file
