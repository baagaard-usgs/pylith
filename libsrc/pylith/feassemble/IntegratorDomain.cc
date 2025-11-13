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

#include "pylith/feassemble/IntegratorDomain.hh" // implementation of object methods

#include "pylith/feassemble/UpdateStateVars.hh" // HOLDSA UpdateStateVars
#include "pylith/feassemble/DSLabelAccess.hh" // USES DSLabelAccess
#include "pylith/problems/Physics.hh" // USES Physics
#include "pylith/feassemble/IntegrationData.hh" // USES IntegrationData
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface::FaceEnum
#include "pylith/feassemble/InterfacePatches.hh" // USES InterfacePatches
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createSubdomainMesh()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField
#include "petscds.h" // USES PetscDS

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error

extern "C" PetscErrorCode DMPlexComputeResidual_Internal(PetscDM dm,
                                                         PetscFormKey key,
                                                         PetscIS cellIS,
                                                         pylith::real time,
                                                         PetscVec locX,
                                                         PetscVec locX_t,
                                                         pylith::real t,
                                                         PetscVec locF,
                                                         void *user);

extern "C" PetscErrorCode DMPlexComputeJacobian_Internal(PetscDM dm,
                                                         PetscFormKey key,
                                                         PetscIS cellIS,
                                                         pylith::real t,
                                                         pylith::real X_tShift,
                                                         PetscVec X,
                                                         PetscVec X_t,
                                                         PetscMat Jac,
                                                         PetscMat JacP,
                                                         void *user);

extern "C" PetscErrorCode DMPlexComputeJacobian_Action_Internal(PetscDM,
                                                                PetscFormKey,
                                                                PetscIS,
                                                                pylith::real,
                                                                pylith::real,
                                                                PetscVec,
                                                                PetscVec,
                                                                PetscVec,
                                                                PetscVec,
                                                                void *);

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace feassemble {
        // Trampoline class for IntegratorDomain::create() factory method
        class IntegratorDomainWrap : public IntegratorDomain {
public:

            IntegratorDomainWrap(std::shared_ptr<pylith::problems::Physics>& physics) : IntegratorDomain(physics) {}


        };

        class _IntegratorDomain {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static pylith::integer initialize;
                static pylith::integer setInterfaceData;
                static pylith::integer setState;
                static pylith::integer computeRHSResidual;
                static pylith::integer computeLHSResidual;
                static pylith::integer computeLHSJacobian;
                static pylith::integer computeLHSJacobianLumpedInv;
                static pylith::integer updateStateVars;
                static pylith::integer computeDerivedField;
            };
        };

    }
}
pylith::utils::EventLogger pylith::feassemble::_IntegratorDomain::Events::logger;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::initialize;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::setInterfaceData;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::setState;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::computeRHSResidual;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::computeLHSResidual;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::computeLHSJacobian;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::computeLHSJacobianLumpedInv;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::updateStateVars;
pylith::integer pylith::feassemble::_IntegratorDomain::Events::computeDerivedField;

// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::_IntegratorDomain::Events::init(void) {
    logger.setClassName("IntegratorDomain");
    logger.initialize();
    initialize = logger.registerEvent("PL:IntegratorDomain:initialize");
    setInterfaceData = logger.registerEvent("PL:IntegratorDomain:setInterfaceData");
    setState = logger.registerEvent("PL:IntegratorDomain:setState");
    computeRHSResidual = logger.registerEvent("PL:IntegratorDomain:computeRHSResidual");
    computeLHSResidual = logger.registerEvent("PL:IntegratorDomain:computeLHSResidual");
    computeLHSJacobian = logger.registerEvent("PL:IntegratorDomain:computeLHSJacobian");
    computeLHSJacobianLumpedInv = logger.registerEvent("PL:IntegratorDomain:computeLHSJacobianLumpedInv");
    updateStateVars = logger.registerEvent("PL:IntegratorDomain:updateStateVars");
    computeDerivedField = logger.registerEvent("PL:IntegratorDomain:computeDerivedField");
}


// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorDomain::IntegratorDomain(std::shared_ptr<pylith::problems::Physics>& physics) :
    Integrator(physics),
    _materialMesh(nullptr),
    _updateState(nullptr),
    _jacobianValues(nullptr),
    _dsLabel(nullptr) {
    GenericComponent::setName("integratordomain");
    _IntegratorDomain::Events::init();
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::feassemble::IntegratorDomain>
pylith::feassemble::IntegratorDomain::create(std::shared_ptr<pylith::problems::Physics>& physics) {
    return std::make_shared<IntegratorDomainWrap>(physics);
} // create


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::feassemble::IntegratorDomain::~IntegratorDomain(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::IntegratorDomain::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::feassemble::Integrator::deallocate();

    _kernelsUpdateStateVars.resize(0);_kernelsUpdateStateVars.shrink_to_fit();
    _kernelsDerivedField.resize(0);_kernelsDerivedField.shrink_to_fit();
    _materialMesh.reset();
    _updateState.reset();
    _jacobianValues.reset();
    _dsLabel.reset();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Get mesh associated with integration domain.
const pylith::topology::Mesh&
pylith::feassemble::IntegratorDomain::getPhysicsDomainMesh(void) const {
    assert(_materialMesh);
    return *_materialMesh;
} // domainMesh


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsResidual(const std::vector<ResidualKernels>& kernels,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsResidual(# kernels="<<kernels.size()<<")");

    PetscErrorCode err = PETSC_SUCCESS;
    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const pylith::integer i_field = solution.getSubfieldInfo(kernels[i].subfield.c_str()).index;
        const pylith::integer i_part = kernels[i].part;
        if (dsLabel.weakForm()) {
            err = PetscWeakFormAddResidual(dsLabel.weakForm(), dsLabel.label(), dsLabel.value(), i_field, i_part,
                                           kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);
        } // if

        switch (kernels[i].part) {
        case LHS:
            _hasLHSResidual = true;
            break;
        case RHS:
            _hasRHSResidual = true;
            break;
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown residual part " << kernels[i].part <<".");
        } // switch
    } // for

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernelsResidual


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsJacobian(const std::vector<JacobianKernels>& kernels,
                                                         const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsJacobian(# kernels="<<kernels.size()<<")");

    PetscErrorCode err = PETSC_SUCCESS;
    DSLabelAccess dsLabel(solution.getDM(), _labelName.c_str(), _labelValue);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const pylith::integer i_fieldTrial = solution.getSubfieldInfo(kernels[i].subfieldTrial.c_str()).index;
        const pylith::integer i_fieldBasis = solution.getSubfieldInfo(kernels[i].subfieldBasis.c_str()).index;
        const pylith::integer i_part = kernels[i].part;
        if (dsLabel.weakForm()) {
            err = PetscWeakFormAddJacobian(dsLabel.weakForm(), dsLabel.label(), dsLabel.value(), i_fieldTrial, i_fieldBasis,
                                           i_part, kernels[i].j0, kernels[i].j1, kernels[i].j2, kernels[i].j3);PYLITH_CHECK_ERROR(err);
        } // if

        switch (kernels[i].part) {
        case LHS:
            _hasLHSJacobian = true;
            break;
        case LHS_LUMPED_INV:
            _hasLHSJacobianLumped = true;
            break;
        default:
            PYLITH_JOURNAL_LOGICERROR("Unknown Jacobian part " << kernels[i].part <<".");
        } // switch
    } // for

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        err = PetscDSView(dsLabel.ds(), PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
} // setKernelsJacobian


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian without finite-element integration.
void
pylith::feassemble::IntegratorDomain::setKernelsJacobian(const std::vector<pylith::feassemble::JacobianValues::JacobianKernel>& kernelsJacobian,
                                                         const std::vector<pylith::feassemble::JacobianValues::JacobianKernel>& kernelsPrecond) {
    _jacobianValues = std::unique_ptr<pylith::feassemble::JacobianValues>();assert(_jacobianValues);
    _jacobianValues->setKernels(kernelsJacobian, kernelsPrecond);
    _hasLHSJacobian = true;
} // setKernelsJacobian


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsUpdateStateVars(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsUpdateStateVars(# kernels="<<kernels.size()<<")");

    _kernelsUpdateStateVars = kernels;

    PYLITH_METHOD_END;
} // setKernelsUpdateStateVars


// ------------------------------------------------------------------------------------------------
void
pylith::feassemble::IntegratorDomain::setKernelsDerivedField(const std::vector<ProjectKernels>& kernels) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setKernelsDerivedField(# kernels="<<kernels.size()<<")");

    _kernelsDerivedField = kernels;

    PYLITH_METHOD_END;
} // setKernelsDerivedField


// ------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorDomain::initialize(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" initialize(solution="<<solution.getName()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::initialize);

    assert(_physics);
    const char* componentName = _physics->getFullIdentifier();
    _materialMesh = pylith::topology::MeshOps::createSubdomainMesh(solution.getMesh(), _labelName.c_str(), _labelValue, componentName);
    pylith::topology::CoordsVisitor::optimizeClosure(_materialMesh->getDM());

    Integrator::initialize(solution);

    assert(_auxiliaryField);
    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmSoln = solution.getDM();assert(dmSoln);
    PetscDMLabel dmLabel = nullptr;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, LHS, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, RHS, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, LHS_LUMPED_INV, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);

    if (_kernelsUpdateStateVars.size() > 0) {
        _updateState = std::unique_ptr<pylith::feassemble::UpdateStateVars>();assert(_updateState);
        _updateState->initialize(*_auxiliaryField);
    } // if

    _dsLabel = std::unique_ptr<DSLabelAccess>(new DSLabelAccess(solution.getDM(), _labelName.c_str(), _labelValue));assert(_dsLabel);
    _dsLabel->removeOverlap();

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing auxiliary field.");
        _auxiliaryField->view("Auxiliary field");
    } // if

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::initialize);
    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Set data needed for integrating faces on interior interfaces.
void
pylith::feassemble::IntegratorDomain::setInterfaceData(const pylith::topology::Field* solution,
                                                       const std::vector<pylith::feassemble::IntegratorInterface*> interfaceIntegrators) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" setInterfaceData(# interfaceIntegrators="<<interfaceIntegrators.size()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::setInterfaceData);
    typedef pylith::feassemble::InterfacePatches::keysmap_t keysmap_t;

    assert(solution);
    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmSoln = solution->getDM();assert(dmSoln);
    PetscDMLabel dmLabel = nullptr;
    err = DMGetLabel(dmSoln, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);

    const pylith::integer numParts = 3;
    const pylith::integer parts[numParts] = {
        LHS,
        RHS,
        LHS_WEIGHTED,
    };
    pylith::integer faultFaces[2];

    const size_t numIntegrators = interfaceIntegrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        const pylith::feassemble::IntegratorInterface* integrator = interfaceIntegrators[i];
        const pylith::feassemble::InterfacePatches& patches = integrator->getIntegrationPatches();
        const keysmap_t patchKeys = patches.getKeys();
        for (keysmap_t::const_iterator iter = patchKeys.begin(); iter != patchKeys.end(); ++iter) {
            size_t faceCount = 0;

            // Check negative face
            const char* negativeLabelName = iter->second.negative.getName();
            const pylith::integer negativeLabelValue = iter->second.negative.getValue();
            if ((std::string(negativeLabelName) == _labelName) && (negativeLabelValue == _labelValue)) {
                faultFaces[faceCount++] = pylith::feassemble::IntegratorInterface::NEGATIVE_FACE;
            } // if

            // Check positive face
            const char* positiveLabelName = iter->second.positive.getName();
            const pylith::integer positiveLabelValue = iter->second.positive.getValue();
            if ((std::string(positiveLabelName) == _labelName) && (positiveLabelValue == _labelValue)) {
                faultFaces[faceCount++] = pylith::feassemble::IntegratorInterface::POSITIVE_FACE;
            } // if

            if (faceCount > 0) { // JOURNAL DEBUGGING
                pythia::journal::debug_t debug(GenericComponent::getName());
                debug << pythia::journal::at(__HERE__) \
                      << "    Found matching interface patch "
                      << patches.getLabelName() << "=" << iter->first << ":";
                for (size_t i = 0; i < faceCount; ++i) {
                    if (faultFaces[i] == pylith::feassemble::IntegratorInterface::NEGATIVE_FACE) {
                        debug << " negative face";
                    } else {
                        if (faultFaces[i] == pylith::feassemble::IntegratorInterface::POSITIVE_FACE) {
                            debug << " positive face";
                        } // if
                    } // if/else
                } // for
                debug << "." << pythia::journal::endl;
            } // JOURNAL DEBUGGING

            const pylith::integer patchValue = iter->second.cohesive.getValue();
            for (size_t iFace = 0; iFace < faceCount; ++iFace) {
                for (pylith::integer iPart = 0; iPart < numParts; ++iPart) {
                    const pylith::integer part = integrator->getWeakFormPart(parts[iPart], faultFaces[iFace], patchValue);

                    assert(_auxiliaryField);
                    err = DMSetAuxiliaryVec(dmSoln, dmLabel, _labelValue, part,
                                            _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
                } // for
            } // for
        } // for
    } // for

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::setInterfaceData);
    PYLITH_METHOD_END;
} // setInterfaceData


// ------------------------------------------------------------------------------------------------
// Set auxiliary field values for current time.
void
pylith::feassemble::IntegratorDomain::setState(const pylith::real t) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setState(t="<<t<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::setState);

    Integrator::setState(t);

    assert(_physics);
    _physics->updateAuxiliaryField(_auxiliaryField.get(), t);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        assert(_auxiliaryField);
        PYLITH_JOURNAL_DEBUG("IntegratorInterface component '" << GenericComponent::getName() << "' for '"
                                                               <<_physics->getIdentifier()
                                                               << "': viewing auxiliary field.");
        _auxiliaryField->view("IntegratorInterface auxiliary field", pylith::topology::Field::VIEW_ALL);
    } // if

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::setState);
    PYLITH_METHOD_END;
} // setState


// ------------------------------------------------------------------------------------------------
// Compute RHS residual for G(t,s).
void
pylith::feassemble::IntegratorDomain::computeRHSResidual(pylith::topology::Field* residual,
                                                         const pylith::feassemble::IntegrationData& integrationData) {
    if (!_hasRHSResidual) { return; }
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeRHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::computeRHSResidual);
    assert(residual);

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::real t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const pylith::real dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);

    _setKernelConstants(*solution, dt);

    assert(_dsLabel);
    PetscFormKey key;
    key.label = _dsLabel->label();
    key.value = _dsLabel->value();
    key.part = pylith::feassemble::Integrator::RHS;

    PetscErrorCode err = PETSC_SUCCESS;
    assert(solution->getLocalVector());
    assert(residual->getLocalVector());
    PetscVec solutionDotVec = nullptr;
    err = DMPlexComputeResidual_Internal(_dsLabel->dm(), key, _dsLabel->cellsIS(), PETSC_MIN_REAL, solution->getLocalVector(),
                                         solutionDotVec, t, residual->getLocalVector(), nullptr);PYLITH_CHECK_ERROR(err);

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::computeRHSResidual);
    PYLITH_METHOD_END;
} // computeRHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS residual for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSResidual(pylith::topology::Field* residual,
                                                         const pylith::feassemble::IntegrationData& integrationData) {
    if (!_hasLHSResidual) { return; }
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSResidual(residual="<<residual<<", integrationData="<<integrationData.str()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::computeLHSResidual);

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::topology::Field* solutionDot = integrationData.getField(pylith::feassemble::IntegrationData::solution_dot);
    assert(solutionDot);
    const pylith::real t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const pylith::real dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);

    _setKernelConstants(*solution, dt);

    assert(_dsLabel);
    PetscFormKey key;
    key.label = _dsLabel->label();
    key.value = _dsLabel->value();
    key.part = pylith::feassemble::Integrator::LHS;

    PetscErrorCode err = PETSC_SUCCESS;
    assert(solution->getLocalVector());
    assert(solutionDot->getLocalVector());
    assert(residual->getLocalVector());
    err = DMPlexComputeResidual_Internal(_dsLabel->dm(), key, _dsLabel->cellsIS(), PETSC_MIN_REAL, solution->getLocalVector(),
                                         solutionDot->getLocalVector(), t, residual->getLocalVector(), nullptr);PYLITH_CHECK_ERROR(err);

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::computeLHSResidual);
    PYLITH_METHOD_END;
} // computeLHSResidual


// ------------------------------------------------------------------------------------------------
// Compute LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSJacobian(PetscMat jacobianMat,
                                                         PetscMat precondMat,
                                                         const pylith::feassemble::IntegrationData& integrationData) {
    if (!_hasLHSJacobian) { return;}
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobian(jacobianMat="<<jacobianMat<<", precondMat="<<precondMat<<", integrationData="<<integrationData.str()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::computeLHSJacobian);

    _needNewLHSJacobian = false;
    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::topology::Field* solutionDot = integrationData.getField(pylith::feassemble::IntegrationData::solution_dot);
    assert(solutionDot);
    const pylith::real t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const pylith::real dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);
    const pylith::real s_tshift = integrationData.getScalar(pylith::feassemble::IntegrationData::s_tshift);

    _setKernelConstants(*solution, dt);

    assert(_dsLabel);
    PetscFormKey key;
    key.label = _dsLabel->label();
    key.value = _dsLabel->value();
    key.part = pylith::feassemble::Integrator::LHS;

    PetscErrorCode err = PETSC_SUCCESS;
    assert(solution->getLocalVector());
    assert(solutionDot->getLocalVector());
    assert(jacobianMat);
    assert(precondMat);
    err = DMPlexComputeJacobian_Internal(_dsLabel->dm(), key, _dsLabel->cellsIS(), t, s_tshift, solution->getLocalVector(),
                                         solutionDot->getLocalVector(), jacobianMat, precondMat, nullptr);PYLITH_CHECK_ERROR(err);

    if (_jacobianValues) {
        _jacobianValues->computeLHSJacobian(jacobianMat, precondMat, t, dt, s_tshift, *solution, *_dsLabel);
        err = MatAssemblyBegin(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
        err = MatAssemblyEnd(jacobianMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
        if (precondMat && ( precondMat != jacobianMat) ) {
            err = MatAssemblyBegin(precondMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
            err = MatAssemblyEnd(precondMat, MAT_FINAL_ASSEMBLY);PYLITH_CHECK_ERROR(err);
        } // if
    } // if

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::computeLHSJacobian);
    PYLITH_METHOD_END;
} // computeLHSJacobian


// ------------------------------------------------------------------------------------------------
// Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}).
void
pylith::feassemble::IntegratorDomain::computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                                                  const pylith::feassemble::IntegrationData& integrationData) {
    if (!_hasLHSJacobianLumped) { return; }
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG(_labelName<<"="<<_labelValue<<" computeLHSJacobianLumpedInv(jacobianInv="<<jacobianInv<<", integrationData="<<integrationData.str()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::computeLHSJacobianLumpedInv);

    _needNewLHSJacobianLumped = false;

    const pylith::topology::Field* solution = integrationData.getField(pylith::feassemble::IntegrationData::solution);
    assert(solution);
    const pylith::real t = integrationData.getScalar(pylith::feassemble::IntegrationData::time);
    const pylith::real dt = integrationData.getScalar(pylith::feassemble::IntegrationData::time_step);
    const pylith::real s_tshift = integrationData.getScalar(pylith::feassemble::IntegrationData::s_tshift);

    _setKernelConstants(*solution, dt);

    assert(_dsLabel);
    PetscFormKey key;
    key.label = _dsLabel->label();
    key.value = _dsLabel->value();
    key.part = pylith::feassemble::Integrator::LHS_LUMPED_INV;

    PetscErrorCode err = PETSC_SUCCESS;
    PetscVec vecRowSum = nullptr;
    err = DMGetLocalVector(_dsLabel->dm(), &vecRowSum);PYLITH_CHECK_ERROR(err);
    err = VecSet(vecRowSum, 1.0);PYLITH_CHECK_ERROR(err);

    assert(jacobianInv);
    assert(jacobianInv->getLocalVector());
    err = DMPlexComputeJacobian_Action_Internal(_dsLabel->dm(), key, _dsLabel->cellsIS(), t, s_tshift, vecRowSum, nullptr,
                                                vecRowSum, jacobianInv->getLocalVector(), nullptr);PYLITH_CHECK_ERROR(err);

    err = DMRestoreLocalVector(_dsLabel->dm(), &vecRowSum);PYLITH_CHECK_ERROR(err);
    // Compute the Jacobian inverse.
    err = VecReciprocal(jacobianInv->getLocalVector());PYLITH_CHECK_ERROR(err);

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::computeLHSJacobianLumpedInv);
    PYLITH_METHOD_END;
} // computeLHSJacobianLumpedInv


// ------------------------------------------------------------------------------------------------
// Update state variables as needed.
void
pylith::feassemble::IntegratorDomain::_updateStateVars(const pylith::real t,
                                                       const pylith::real dt,
                                                       const pylith::topology::Field& solution) {
    if (0 == _kernelsUpdateStateVars.size()) { return; }
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_updateStateVars(t="<<t<<", dt="<<dt<<", solution="<<solution.getName()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::updateStateVars);

    assert(_updateState);
    assert(_auxiliaryField);
    _updateState->prepare(_auxiliaryField.get());
    _setKernelConstants(solution, dt);

    // We assume order of the update state variable kernels matches
    // the order of the correspoinding subfields in the auxiliary
    // field.
    const size_t numKernels = _kernelsUpdateStateVars.size();
    PetscPointFunc* kernelsStateVars = (numKernels > 0) ? new PetscPointFunc[numKernels] : nullptr;
    for (size_t iKernel = 0; iKernel < numKernels; ++iKernel) {
        kernelsStateVars[iKernel] = _kernelsUpdateStateVars[iKernel].f;
    } // for

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM stateVarsDM = _updateState->stateVarsDM();
    PetscDMLabel dmLabel = nullptr;
    pylith::integer labelValue = 0;
    const pylith::integer part = 0;
    err = DMSetAuxiliaryVec(stateVarsDM, dmLabel, labelValue, part, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(stateVarsDM, t, solution.getLocalVector(), kernelsStateVars, INSERT_VALUES,
                              _updateState->stateVarsLocalVector());PYLITH_CHECK_ERROR(err);
    _updateState->restore(_auxiliaryField.get());

    delete[] kernelsStateVars;kernelsStateVars = nullptr;

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::updateStateVars);
    PYLITH_METHOD_END;
} // _updateStateVars


// ------------------------------------------------------------------------------------------------
// Compute field derived from solution and auxiliary field.
void
pylith::feassemble::IntegratorDomain::_computeDerivedField(const pylith::real t,
                                                           const pylith::real dt,
                                                           const pylith::topology::Field& solution) {
    if (!_derivedField) { return; }
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeDerivedField(t="<<t<<", dt="<<dt<<", solution="<<solution.getName()<<")");
    _IntegratorDomain::Events::logger.eventBegin(_IntegratorDomain::Events::computeDerivedField);

    assert(_derivedField);
    _setKernelConstants(solution, dt);

    const size_t numKernels = _kernelsDerivedField.size();
    std::vector<PetscPointFunc> kernelsArray(numKernels);
    for (size_t iKernel = 0; iKernel < numKernels; ++iKernel) {
        const pylith::topology::Field::SubfieldInfo& sinfo = _derivedField->getSubfieldInfo(_kernelsDerivedField[iKernel].subfield.c_str());
        kernelsArray[sinfo.index] = _kernelsDerivedField[iKernel].f;
    } // for

    PetscErrorCode err = PETSC_SUCCESS;

    PetscDM derivedDM = _derivedField->getDM();
    assert(_auxiliaryField);
    PetscDMLabel dmLabel = nullptr;
    pylith::integer labelValue = 0;
    const pylith::integer part = 0;
    err = DMSetAuxiliaryVec(derivedDM, dmLabel, labelValue, part, _auxiliaryField->getLocalVector());PYLITH_CHECK_ERROR(err);
    err = DMProjectFieldLocal(derivedDM, t, solution.getLocalVector(), kernelsArray.data(), INSERT_VALUES, _derivedField->getLocalVector());PYLITH_CHECK_ERROR(err);

    pythia::journal::debug_t debug(GenericComponent::getName());
    if (debug.state()) {
        PYLITH_JOURNAL_DEBUG("Viewing derived field.");
        _derivedField->view("Derived field");
    } // if

    _IntegratorDomain::Events::logger.eventEnd(_IntegratorDomain::Events::computeDerivedField);
    PYLITH_METHOD_END;
} // _computeDerivedField


// End of file
