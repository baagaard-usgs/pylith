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

#include "pylith/problems/Problem.hh" // implementation of class methods

#include "pylith/feassemble/IntegrationData.hh" // HOLDSA IntegrationData
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // HASA Field
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/materials/Material.hh" // USES Material
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/bc/BoundaryCondition.hh" // USES BoundaryCondition
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/feassemble/IntegratorInterface.hh" // USES IntegratorInterface
#include "pylith/feassemble/Constraint.hh" // USES Constraint
#include "pylith/problems/ObserversSoln.hh" // USES ObserversSoln
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/EventLogger.hh" // HASA EventLogger
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <cassert> // USES assert()
#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace problems {
        class _Problem {
public:

            /** Create null space for solution subfield.
             *
             * @param[inout] solution Solution field.
             * @param[in] subfieldName Name of solution subfield with null space.
             */
            static
            void createNullSpace(const pylith::topology::Field* solution,
                                 const char* subfieldName);

            /** Set data needed to integrate domain faces on interior interface.
             *
             * @param[inout] solution Solution field.
             * @param[in] integrators Array of integrators for problem.
             */
            static
            void setInterfaceData(const pylith::topology::Field* solution,
                                  std::vector<pylith::feassemble::Integrator*>& integrators);

            /** Get subset of integrators matching template type T.
             *
             * @param[in] Array of integrators for problem
             * @returns Array of integrators of type T.
             */
            template<class T> static std::vector<T*> subset(const std::vector<pylith::feassemble::Integrator*>& integrators);

            // Logging events
            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static pylith::integer setSolution;
                static pylith::integer preinitialize;
                static pylith::integer verifyConfiguration;
                static pylith::integer initialize;
                static pylith::integer checkMaterials;
                static pylith::integer createIntegrators;
                static pylith::integer createConstraints;
                static pylith::integer setupSolution;
            };

        };
    }
}
pylith::utils::EventLogger pylith::problems::_Problem::Events::logger;
pylith::integer pylith::problems::_Problem::Events::setSolution;
pylith::integer pylith::problems::_Problem::Events::preinitialize;
pylith::integer pylith::problems::_Problem::Events::verifyConfiguration;
pylith::integer pylith::problems::_Problem::Events::initialize;
pylith::integer pylith::problems::_Problem::Events::checkMaterials;
pylith::integer pylith::problems::_Problem::Events::createIntegrators;
pylith::integer pylith::problems::_Problem::Events::createConstraints;
pylith::integer pylith::problems::_Problem::Events::setupSolution;

// ------------------------------------------------------------------------------------------------
void
pylith::problems::_Problem::Events::init(void) {
    logger.setClassName("Problem");
    logger.initialize();
    setSolution = logger.registerEvent("PL:Problem:setSolution");
    preinitialize = logger.registerEvent("PL:Problem:preinitialize");
    verifyConfiguration = logger.registerEvent("PL:Problem:verifyConfiguration");
    initialize = logger.registerEvent("PL:Problem:initialize");
    checkMaterials = logger.registerEvent("PL:Problem:initialize");
    createIntegrators = logger.registerEvent("PL:Problem:initialize");
    createConstraints = logger.registerEvent("PL:Problem:initialize");
    setupSolution = logger.registerEvent("PL:Problem:initialize");
} // init


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::problems::Problem::Problem() :
    _integrationData(new pylith::feassemble::IntegrationData),
    _normalizer(nullptr),
    _gravityField(nullptr),
    _observers(new pylith::problems::ObserversSoln),
    _formulation(pylith::problems::Physics::QUASISTATIC),
    _solverType(LINEAR),
    _petscDefaults(pylith::utils::PetscDefaults::SOLVER | pylith::utils::PetscDefaults::TESTING) {
    _Problem::Events::init();
}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::Problem::~Problem(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::Problem::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _integrationData.reset();
    _normalizer.reset();
    _gravityField.reset();
    _materials.resize(0);_materials.shrink_to_fit();
    _bcs.resize(0);_bcs.shrink_to_fit();
    _interfaces.resize(0);_interfaces.shrink_to_fit();
    _integrators.resize(0);_integrators.shrink_to_fit();
    _constraints.resize(0);_constraints.shrink_to_fit();
    _observers.reset();

    pylith::topology::FieldOps::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Set formulation for solving equation.
void
pylith::problems::Problem::setFormulation(const pylith::problems::Physics::FormulationEnum value) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setFormulation(value="<<value<<")");

    _formulation = value;

    PYLITH_METHOD_END;
} // setFormulation


// ------------------------------------------------------------------------------------------------
// Get formulation for solving equation.
pylith::problems::Physics::FormulationEnum
pylith::problems::Problem::getFormulation(void) const {
    return _formulation;
} // getFormulation


// ------------------------------------------------------------------------------------------------
// Set problem type.
void
pylith::problems::Problem::setSolverType(const SolverTypeEnum value) {
    PYLITH_COMPONENT_DEBUG("Problem::setSolverType(value="<<value<<")");

    _solverType = value;
} // setSolverType


// ------------------------------------------------------------------------------------------------
// Get problem type.
pylith::problems::Problem::SolverTypeEnum
pylith::problems::Problem::getSolverType(void) const {
    return _solverType;
} // getSolverType


// ------------------------------------------------------------------------------------------------
// Specify whether to set defaults for PETSc solver appropriate for problem.
void
pylith::problems::Problem::setPetscDefaults(const int flags) {
    _petscDefaults = flags;
} // setPetscDefaults


// ------------------------------------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::problems::Problem::setNormalizer(const std::shared_ptr<spatialdata::units::Nondimensional>& normalizer) {
    PYLITH_COMPONENT_DEBUG("Problem::setNormalizer(dim="<<typeid(normalizer).name()<<")");

    _normalizer = normalizer;
} // setNormalizer


// ------------------------------------------------------------------------------------------------
// Set gravity field.
void
pylith::problems::Problem::setGravityField(const std::shared_ptr<spatialdata::spatialdb::GravityField>& gravityField) {
    PYLITH_COMPONENT_DEBUG("Problem::setGravityField(g="<<typeid(gravityField).name()<<")");

    _gravityField = gravityField;
} // setGravityField


// ----------------------------------------------------------------------
// Register observer to receive notifications.
void
pylith::problems::Problem::registerObserver(const std::shared_ptr<pylith::problems::ObserverSoln>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("registerObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    assert(_normalizer);
    _observers->registerObserver(observer);
    _observers->setTimeScale(_normalizer->getTimeScale());

    PYLITH_METHOD_END;
} // registerObserver


// ----------------------------------------------------------------------
// Remove observer from receiving notifications.
void
pylith::problems::Problem::removeObserver(const std::shared_ptr<pylith::problems::ObserverSoln>& observer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("removeObserver(observer="<<typeid(observer).name()<<")");

    assert(_observers);
    _observers->removeObserver(observer);

    PYLITH_METHOD_END;
} // removeObserver


// ------------------------------------------------------------------------------------------------
// Set solution field.
void
pylith::problems::Problem::setSolution(const std::shared_ptr<pylith::topology::Field>& field) {
    PYLITH_COMPONENT_DEBUG("Problem::setSolution(field="<<typeid(*field).name()<<")");
    _Problem::Events::logger.eventBegin(_Problem::Events::setSolution);

    assert(_integrationData);
    _integrationData->setField(pylith::feassemble::IntegrationData::solution, field);

    _Problem::Events::logger.eventEnd(_Problem::Events::setSolution);
} // setSolution


// ------------------------------------------------------------------------------------------------
// Get solution field.
const pylith::topology::Field*
pylith::problems::Problem::getSolution(void) const {
    PYLITH_METHOD_BEGIN;

    assert(_integrationData);
    pylith::topology::Field* solution = nullptr;
    if (_integrationData->hasField(pylith::feassemble::IntegrationData::solution)) {
        solution = _integrationData->getField(pylith::feassemble::IntegrationData::solution);
    } // if

    PYLITH_METHOD_RETURN(solution);
} // getSolution


// ------------------------------------------------------------------------------------------------
// Get time derivative solution field.
const pylith::topology::Field*
pylith::problems::Problem::getSolutionDot(void) const {
    assert(_integrationData);
    pylith::topology::Field* solutionDot = nullptr;
    if (_integrationData->hasField(pylith::feassemble::IntegrationData::solution_dot)) {
        solutionDot = _integrationData->getField(pylith::feassemble::IntegrationData::solution_dot);
    } // if

    PYLITH_METHOD_RETURN(solutionDot);
} // getSolutionDot


// ------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setMaterials(const std::vector<std::shared_ptr<pylith::materials::Material> >& materials) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setMaterials("<<materials<<")");

    _materials = materials;

    PYLITH_METHOD_END;
} // setMaterials


// ------------------------------------------------------------------------------------------------
// Set boundary conditions.
void
pylith::problems::Problem::setBoundaryConditions(const std::vector<std::shared_ptr<pylith::bc::BoundaryCondition> >& bcs) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setBoundaryConditions("<<bcs<<")");

    _bcs = bcs;

    PYLITH_METHOD_END;
} // setBoundaryConditions


// ------------------------------------------------------------------------------------------------
// Set materials.
void
pylith::problems::Problem::setInterfaces(const std::vector<std::shared_ptr<pylith::faults::FaultCohesive> >& interfaces) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::setInterfaces("<<interfaces<<")");

    _interfaces = interfaces;

    PYLITH_METHOD_END;
} // setInterfaces


// ----------------------------------------------------------------------
// Do minimal initialization.
void
pylith::problems::Problem::preinitialize(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::preinitialzie(mesh="<<typeid(mesh).name()<<")");
    _Problem::Events::logger.eventBegin(_Problem::Events::preinitialize);

    assert(_normalizer);

    for (auto material : _materials) {
        assert(material);
        material->setNormalizer(_normalizer);
        material->setGravityField(_gravityField);
        material->setFormulation(_formulation);
    } // for

    for (auto interface : _interfaces) {
        assert(interface);
        interface->setNormalizer(_normalizer);
        interface->setFormulation(_formulation);
    } // for

    for (auto bc : _bcs) {
        assert(bc);
        bc->setNormalizer(_normalizer);
        bc->setFormulation(_formulation);
    } // for

    _Problem::Events::logger.eventEnd(_Problem::Events::preinitialize);
    PYLITH_METHOD_END;
} // preinitialize


// ----------------------------------------------------------------------
// Verify configuration.
void
pylith::problems::Problem::verifyConfiguration(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::verifyConfiguration(void)");
    _Problem::Events::logger.eventBegin(_Problem::Events::verifyConfiguration);

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    // Check to make sure materials are compatible with the solution.
    for (auto material : _materials) {
        assert(material);
        material->verifyConfiguration(*solution);
    } // for

    // Check to make sure interfaces are compatible with the solution.
    for (auto interface : _interfaces) {
        assert(interface);
        interface->verifyConfiguration(*solution);
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    for (auto bc : _bcs) {
        assert(bc);
        bc->verifyConfiguration(*solution);
    } // for

    _checkMaterialLabels();

    assert(_observers);
    _observers->verifyObservers(*solution);

    _Problem::Events::logger.eventEnd(_Problem::Events::verifyConfiguration);
    PYLITH_METHOD_END;
} // verifyConfiguration


// ----------------------------------------------------------------------
// Initialize.
void
pylith::problems::Problem::initialize(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::initialize()");
    _Problem::Events::logger.eventBegin(_Problem::Events::initialize);

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    // Initialize solution field.
    pylith::utils::PetscDefaults::set(*solution, _materials[0], _petscDefaults);
    PetscErrorCode err = DMSetFromOptions(solution->getDM());PYLITH_CHECK_ERROR(err);
    _setupSolution();
    pylith::topology::CoordsVisitor::optimizeClosure(solution->getDM());

    // Initialize integrators.
    _createIntegrators();
    const size_t numIntegrators = _integrators.size();
    for (size_t i = 0; i < numIntegrators; ++i) {
        assert(_integrators[i]);
        _integrators[i]->initialize(*solution);
    } // for

    // Initialize constraints.
    _createConstraints();
    const size_t numConstraints = _constraints.size();
    for (size_t i = 0; i < numConstraints; ++i) {
        assert(_constraints[i]);
        _constraints[i]->initialize(*solution);
    } // for

    solution->allocate();
    solution->createGlobalVector();
    solution->createOutputVector();

    switch (_formulation) {
    case pylith::problems::Physics::DYNAMIC:
    case pylith::problems::Physics::DYNAMIC_IMEX:
        break;
    case pylith::problems::Physics::QUASISTATIC:
        _Problem::createNullSpace(solution, "displacement");
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '"<<_formulation<<".");
    } // switch
    _Problem::setInterfaceData(solution, _integrators);

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying solution field layout");
        solution->view("Solution field", pylith::topology::Field::VIEW_LAYOUT);
    } // if

    _Problem::Events::logger.eventEnd(_Problem::Events::initialize);
    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Check material and interface ids.
void
pylith::problems::Problem::_checkMaterialLabels(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_checkMaterialLabels()");
    _Problem::Events::logger.eventBegin(_Problem::Events::checkMaterials);

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();

    pylith::integer_array labelValues(numMaterials + numInterfaces);
    size_t count = 0;
    for (size_t i = 0; i < numMaterials; ++i) {
        assert(_materials[i]);
        labelValues[count++] = _materials[i]->getLabelValue();
    } // for
    for (size_t i = 0; i < numInterfaces; ++i) {
        assert(_interfaces[i]);
        labelValues[count++] = _interfaces[i]->getCohesiveLabelValue();
    } // for

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);
    pylith::topology::MeshOps::checkMaterialLabels(solution->getMesh(), labelValues);

    _Problem::Events::logger.eventEnd(_Problem::Events::checkMaterials);
    PYLITH_METHOD_END;
} // _checkMaterialLabels


// ------------------------------------------------------------------------------------------------
// Create array of integrators from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createIntegrators(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_createIntegrators()");
    _Problem::Events::logger.eventBegin(_Problem::Events::createIntegrators);

    const size_t numMaterials = _materials.size();
    const size_t numInterfaces = _interfaces.size();
    const size_t numBC = _bcs.size();
    const size_t maxSize = numMaterials + numInterfaces + numBC;
    _integrators.resize(maxSize);
    size_t count = 0;

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    for (auto material : _materials) {
        assert(material);
        std::unique_ptr<pylith::feassemble::Integrator> integrator = material->createIntegrator(*solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = std::move(integrator);}
    } // for

    for (auto interface : _interfaces) {
        assert(interface);
        std::unique_ptr<pylith::feassemble::Integrator> integrator = interface->createIntegrator(*solution, _materials);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = std::move(integrator);}
    } // for

    // Check to make sure boundary conditions are compatible with the solution.
    for (auto bc : _bcs) {
        assert(bc);
        pylith::feassemble::Integrator* integrator = bc->createIntegrator(*solution);
        assert(count < maxSize);
        if (integrator) { _integrators[count++] = std::move(integrator);}
    } // for

    _integrators.resize(count);

    _Problem::Events::logger.eventEnd(_Problem::Events::createIntegrators);
    PYLITH_METHOD_END;
} // _createIntegrators


// ------------------------------------------------------------------------------------------------
// Create array of constraints from materials, interfaces, and boundary conditions.
void
pylith::problems::Problem::_createConstraints(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_createConstraints()");
    _Problem::Events::logger.eventBegin(_Problem::Events::createConstraints);

    assert(_integrationData);
    const pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);

    _constraints.resize(0); // insure we start with an empty array.

    for (auto material : _materials) {
        assert(material);
        std::vector<pylith::feassemble::Constraint*> constraints = material->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());
    } // for

    for (auto interface : _interfaces) {
        assert(interface);
        std::vector<pylith::feassemble::Constraint*> constraints = interface->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());
    } // for

    for (auto bc : _bcs) {
        assert(bc);
        std::vector<pylith::feassemble::Constraint*> constraints = bc->createConstraints(*solution);
        _constraints.insert(_constraints.end(), constraints.begin(), constraints.end());
    } // for

    _Problem::Events::logger.eventEnd(_Problem::Events::createConstraints);
    PYLITH_METHOD_END;
} // _createConstraints


// ------------------------------------------------------------------------------------------------
// Setup solution subfields and discretization.
void
pylith::problems::Problem::_setupSolution(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("Problem::_setupSolution()");
    _Problem::Events::logger.eventBegin(_Problem::Events::createConstraints);

    assert(_integrationData);
    pylith::topology::Field* solution = _integrationData->getField("solution");
    assert(solution);
    solution->subfieldsSetup();
    solution->createDiscretization();

    // Mark fault fields as implicit.
    const pylith::string_vector& subfieldNames = solution->getSubfieldNames();
    for (size_t i = 0; i < subfieldNames.size(); ++i) {
        const pylith::topology::Field::SubfieldInfo& subfieldInfo = solution->getSubfieldInfo(subfieldNames[i].c_str());
        if (subfieldInfo.fe.isFaultOnly) {
            PetscErrorCode err = PETSC_SUCCESS;
            PetscDS ds = nullptr;
            pylith::integer cStart = 0, cEnd = 0;
            PetscDM dmSoln = solution->getDM();assert(dmSoln);
            err = DMPlexGetHeightStratum(dmSoln, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
            pylith::integer cell = cStart;
            bool found = false;
            for (; cell < cEnd; ++cell) {
                if (pylith::topology::MeshOps::isCohesiveCell(dmSoln, cell)) {
                    found = true;
                    break;
                } // if
            } // for
            if (!found) {
                continue;
            } // if
            err = DMGetCellDS(dmSoln, cell, &ds, nullptr);PYLITH_CHECK_ERROR(err);
            assert(ds);
            err = PetscDSSetImplicit(ds, subfieldInfo.index, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    _Problem::Events::logger.eventEnd(_Problem::Events::setupSolution);
    PYLITH_METHOD_END;
} // _setupSolution


// ------------------------------------------------------------------------------------------------
// Create null space for solution subfield.
void
pylith::problems::_Problem::createNullSpace(const pylith::topology::Field* solution,
                                            const char* subfieldName) {
    PYLITH_METHOD_BEGIN;
    assert(solution);

    const int spaceDim = solution->getSpaceDim();
    const pylith::integer m = (spaceDim * (spaceDim + 1)) / 2;assert(m > 0 && m <= 6);

    PetscErrorCode err = PETSC_SUCCESS;
    pylith::integer numDofUnconstrained = 0;
    err = PetscSectionGetConstrainedStorageSize(solution->getLocalSection(), &numDofUnconstrained);
    if (m > numDofUnconstrained) {
        PYLITH_METHOD_END;
    } // if

    const PetscDM dmSoln = solution->getDM();
    const pylith::topology::Field::SubfieldInfo info = solution->getSubfieldInfo(subfieldName);
    MatNullSpace nullSpace = nullptr;
    err = DMPlexCreateRigidBody(dmSoln, info.index, &nullSpace);PYLITH_CHECK_ERROR(err);

    PetscObject field = nullptr;
    err = DMGetField(dmSoln, info.index, nullptr, &field);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose(field, "nearnullspace", (PetscObject) nullSpace);PYLITH_CHECK_ERROR(err);
    err = MatNullSpaceDestroy(&nullSpace);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // createNullSpace


// ------------------------------------------------------------------------------------------------
// Set data needed to integrate domain faces on interior interface.
void
pylith::problems::_Problem::setInterfaceData(const pylith::topology::Field* solution,
                                             std::vector<pylith::feassemble::Integrator*>& integrators) {
    PYLITH_METHOD_BEGIN;

    const std::vector<pylith::feassemble::IntegratorDomain*>& integratorsDomain =
        subset<pylith::feassemble::IntegratorDomain>(integrators);
    const std::vector<pylith::feassemble::IntegratorInterface*>& integratorsInterface =
        subset<pylith::feassemble::IntegratorInterface>(integrators);

    for (size_t i = 0; i < integratorsDomain.size(); ++i) {
        integratorsDomain[i]->setInterfaceData(solution, integratorsInterface);
    } // for

    PYLITH_METHOD_END;
} // setInterfaceData


// ------------------------------------------------------------------------------------------------
// Get subset of integrators matching template type T.
template<class T>
std::vector<T*>
pylith::problems::_Problem::subset(const std::vector<pylith::feassemble::Integrator*>& integrators) {
    // Count number of matches
    size_t numFound = 0;
    for (size_t i = 0; i < integrators.size(); ++i) {
        if (dynamic_cast<T*>(integrators[i])) {
            numFound++;
        } // if/else
    } // for

    // Collect matches
    std::vector<T*> matches(numFound);
    size_t index = 0;
    for (size_t i = 0; i < integrators.size(); ++i) {
        T* integrator = dynamic_cast<T*>(integrators[i]);
        if (integrator) {
            matches[index++] = integrator;
        } // if/else
    } // for

    return matches;
} // subset


// End of file
