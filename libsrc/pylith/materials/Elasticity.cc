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

#include "pylith/materials/Elasticity.hh" // implementation of object methods

#include "pylith/materials/RheologyElasticity.hh" // HASA RheologyElasticity
#include "pylith/materials/SubfieldFactory.hh" // USES SubfieldFactory
#include "pylith/feassemble/IntegratorDomain.hh" // USES IntegratorDomain
#include "pylith/feassemble/JacobianValues.hh" // USES JacobianValues
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field::SubfieldInfo
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/FieldOps.hh" // USES FieldOps

#include "pylith/utils/PetscOptions.hh" // USES PetscOptions

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels
#include "pylith/fekernels/DispVel.hh" // USES DispVel kernels
#include "pylith/fekernels/FaultCohesiveKin.hh" // USES FaultCohesiveKin kernels

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/EventLogger.hh" // EventLogger

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace materials {
        class ElasticityWrap : public Elasticity {
public:

            ElasticityWrap(void) : Elasticity() {}


        };

        class _Elasticity {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static pylith::integer verifyConfiguration;
                static pylith::integer createIntegrator;
                static pylith::integer createAuxiliaryField;
                static pylith::integer createDerivedField;
            };

        }; // _Elasticity
    } // materials
} // pylith

pylith::utils::EventLogger pylith::materials::_Elasticity::Events::logger;
pylith::integer pylith::materials::_Elasticity::Events::verifyConfiguration;
pylith::integer pylith::materials::_Elasticity::Events::createIntegrator;
pylith::integer pylith::materials::_Elasticity::Events::createAuxiliaryField;
pylith::integer pylith::materials::_Elasticity::Events::createDerivedField;

// ------------------------------------------------------------------------------------------------
void
pylith::materials::_Elasticity::Events::init(void) {
    logger.setClassName("Elasticity");
    logger.initialize();
    verifyConfiguration = logger.registerEvent("PL:Elasticity:verifyConfiguration");
    createIntegrator = logger.registerEvent("PL:Elasticity:createIntegrator");
    createAuxiliaryField = logger.registerEvent("PL:Elasticity:createAuxiliaryField");
    createDerivedField = logger.registerEvent("PL:Elasticity:createDerivedField");
}


// ------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ResidualKernels ResidualKernels;
typedef pylith::feassemble::IntegratorDomain::JacobianKernels JacobianKernels;
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;
typedef pylith::feassemble::Integrator::EquationPart EquationPart;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::Elasticity::Elasticity(void) :
    _useBodyForce(false),
    _rheology(nullptr) {
    pylith::utils::PyreComponent::setName("elasticity");
    _Elasticity::Events::init();
} // constructor


// ------------------------------------------------------------------------------------------------
// Factory for std::shared_ptr.
std::shared_ptr<pylith::materials::Elasticity>
pylith::materials::Elasticity::create(void) {
    return std::make_shared<ElasticityWrap>();
}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::Elasticity::~Elasticity(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::Elasticity::deallocate(void) {
    Material::deallocate();

    _rheology.reset();
} // deallocate


// ------------------------------------------------------------------------------------------------
// Include body force?
void
pylith::materials::Elasticity::useBodyForce(const bool value) {
    PYLITH_COMPONENT_DEBUG("useBodyForce(value="<<value<<")");

    _useBodyForce = value;
} // useBodyForce


// ------------------------------------------------------------------------------------------------
// Include body force?
bool
pylith::materials::Elasticity::useBodyForce(void) const {
    return _useBodyForce;
} // useBodyForce


// ------------------------------------------------------------------------------------------------
// Set bulk rheology.
void
pylith::materials::Elasticity::setBulkRheology(const std::shared_ptr<pylith::materials::RheologyElasticity>& rheology) {
    _rheology = rheology;
} // setBulkRheology


// ------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::materials::Elasticity::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getName()<<")");
    _Elasticity::Events::logger.eventBegin(_Elasticity::Events::verifyConfiguration);

    // Verify solution contains required fields.
    std::string reason = "material 'Elasticity'.";
    size_t numRequired = 0;
    const size_t maxRequired = 2;
    pylith::string_vector requiredFields(maxRequired);
    requiredFields[numRequired++] = "displacement";

    switch (_formulation) {
    case QUASISTATIC:
        break;
    case DYNAMIC:
    case DYNAMIC_IMEX:
        reason = "material 'Elasticity' with inertia";
        requiredFields[numRequired++] = "velocity";
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch
    requiredFields.resize(numRequired);

    pylith::topology::FieldOps::checkSubfieldsExist(requiredFields, reason, solution);

    _Elasticity::Events::logger.eventEnd(_Elasticity::Events::verifyConfiguration);
    PYLITH_METHOD_END;
} // verifyConfiguration


// ------------------------------------------------------------------------------------------------
// Create integrator and set kernels.
std::shared_ptr<pylith::feassemble::Integrator>
pylith::materials::Elasticity::createIntegrator(const pylith::topology::Field& solution) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createIntegrator(solution="<<solution.getName()<<")");
    _Elasticity::Events::logger.eventBegin(_Elasticity::Events::createIntegrator);

    std::shared_ptr<pylith::feassemble::IntegratorDomain> integrator = pylith::feassemble::IntegratorDomain::create(shared_from_this());assert(integrator);
    integrator->setLabelName(getLabelName());
    integrator->setLabelValue(getLabelValue());

    _setKernelsResidual(integrator.get(), solution);
    _setKernelsJacobian(integrator.get(), solution);
    _setKernelsUpdateStateVars(integrator.get(), solution);
    _setKernelsDerivedField(integrator.get(), solution);

    _Elasticity::Events::logger.eventEnd(_Elasticity::Events::createIntegrator);
    PYLITH_METHOD_RETURN(integrator);
} // createIntegrator


// ------------------------------------------------------------------------------------------------
// Create auxiliary field.
std::shared_ptr<pylith::topology::Field>
pylith::materials::Elasticity::createAuxiliaryField(const pylith::topology::Field& solution,
                                                    const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createAuxiliaryField(solution="<<solution.getName()<<", domainMesh="<<typeid(domainMesh).name()<<")");
    _Elasticity::Events::logger.eventBegin(_Elasticity::Events::createAuxiliaryField);

    std::shared_ptr<pylith::topology::Field> auxiliaryField = std::make_shared<pylith::topology::Field>(domainMesh);assert(auxiliaryField);
    auxiliaryField->setName("auxiliary field");

    assert(_subfieldFactory);
    assert(_normalizer);
    assert(_rheology);

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the FE kernels.
    //
    // :ATTENTION: In quasi-static problems, the time scale is usually quite large
    // (order of tens to hundreds of years), which means that the density scale is very large,
    // and the acceleration scale is very small. Nevertheless, density times gravitational
    // acceleration will have a scale of pressure divided by length and should be within a few orders
    // of magnitude of 1.
    _subfieldFactory->open(auxiliaryField, *_normalizer);
    _subfieldFactory->addSubfield(_subfieldFactory->density);
    if (_useBodyForce) {
        _subfieldFactory->addSubfield(_subfieldFactory->body_force);
    } // if
    if (_gravityField) {
        _subfieldFactory->addSubfield(_subfieldFactory->gravity_field);
    } // if
    _rheology->addAuxiliarySubfields(_subfieldFactory);
    _subfieldFactory->close();

    auxiliaryField->subfieldsSetup();
    auxiliaryField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *auxiliaryField);
    auxiliaryField->allocate();
    auxiliaryField->createOutputVector();

    pylith::topology::FieldQuery fieldQuery(auxiliaryField);
    fieldQuery.addSubfield(_subfieldFactory->density);
    if (_useBodyForce) {
        fieldQuery.addSubfield(_subfieldFactory->body_force);
    } // if
    if (_gravityField) {
        std::shared_ptr<spatialdata::spatialdb::SpatialDB> db = _gravityField;
        fieldQuery.addSubfield(_subfieldFactory->gravity_field, db);
    } // if
    _rheology->addAuxiliarySubfields(fieldQuery);
    fieldQuery.open(_auxiliaryFieldDB, _normalizer->getLengthScale());
    fieldQuery.query();
    fieldQuery.close();

    _Elasticity::Events::logger.eventEnd(_Elasticity::Events::createAuxiliaryField);
    PYLITH_METHOD_RETURN(auxiliaryField);
} // createAuxiliaryField


// ------------------------------------------------------------------------------------------------
// Create derived field.
std::shared_ptr<pylith::topology::Field>
pylith::materials::Elasticity::createDerivedField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& domainMesh) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("createDerivedField(solution="<<solution.getName()<<", domainMesh="<<typeid(domainMesh).name()<<")");
    _Elasticity::Events::logger.eventBegin(_Elasticity::Events::createDerivedField);

    assert(_subfieldFactory);

    std::shared_ptr<pylith::topology::Field> derivedField = std::make_shared<pylith::topology::Field>(domainMesh);assert(derivedField);
    derivedField->setName("derived field");

    assert(_normalizer);
    _subfieldFactory->open(derivedField, _normalizer);
    pylith::string_vector subfields = {
        _subfieldFactory->cauchy_stress,
        _subfieldFactory->cauchy_strain,
    };
    size_t numSubfields = 0;
    for (auto subfieldName : subfields) {
        if (_subfieldFactory->hasDiscretization(subfieldName.c_str())) {
            _subfieldFactory->addSubfield(subfieldName.c_str());
            ++numSubfields;
        } // if
    } // for
    _subfieldFactory->close();
    if (!numSubfields) {
        derivedField.reset();
    } // if

    derivedField->subfieldsSetup();
    derivedField->createDiscretization();
    pylith::topology::FieldOps::checkDiscretization(solution, *derivedField);
    derivedField->allocate();
    derivedField->createOutputVector();

    _Elasticity::Events::logger.eventEnd(_Elasticity::Events::createDerivedField);
    PYLITH_METHOD_RETURN(derivedField);
} // createDerivedField


// ------------------------------------------------------------------------------------------------
// Get default PETSc solver options appropriate for material.
std::unique_ptr<pylith::utils::PetscOptions>
pylith::materials::Elasticity::getSolverDefaults(const bool isParallel,
                                                 const bool hasFault) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getSolverDefaults(isParallel="<<isParallel<<", hasFault="<<hasFault<<")");

    std::unique_ptr<pylith::utils::PetscOptions> options(new pylith::utils::PetscOptions);assert(options);

    switch (_formulation) {
    case pylith::problems::Physics::QUASISTATIC:
        options->add("-ts_type", "beuler");

        if (!hasFault) {
            if (!isParallel) {
                options->add("-pc_type", "lu");
            } else {
                options->add("-pc_type", "gamg");
            } // if/else
        } else {
            options->add("-pc_type", "gamg");
            options->add("-dm_reorder_section");
            options->add("-dm_reorder_section_type", "cohesive");
            options->add("-mg_fine_pc_type", "vpbjacobi");
        } // if/else
        break;
    case pylith::problems::Physics::DYNAMIC:
        break;
    case pylith::problems::Physics::DYNAMIC_IMEX:
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation '" << _formulation << "'.");
    } // switch

    PYLITH_METHOD_RETURN(options);
} // getSolverDefaults


// ------------------------------------------------------------------------------------------------
// Get residual kernels for an interior interface bounding material.
std::vector<pylith::materials::Material::InterfaceResidualKernels>
pylith::materials::Elasticity::getInterfaceKernelsResidual(const pylith::topology::Field& solution,
                                                           const pylith::feassemble::IntegratorInterface::FaceEnum face) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelsResidual(solution="<<solution.getName()<<", face="<<face<<")");

    const std::shared_ptr<spatialdata::geocoords::CoordSys>& coordsys = solution.getMesh().getCoordSys();

    std::vector<InterfaceResidualKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC:
    case DYNAMIC:
        break;
    case DYNAMIC_IMEX: {
        PetscBdPointFunc f0l = nullptr;
        PetscBdPointFunc f1l = nullptr;

        switch (face) {
        case pylith::feassemble::IntegratorInterface::NEGATIVE_FACE:
            f0l = _rheology->getKernelf0Neg(coordsys.get());
            break;
        case pylith::feassemble::IntegratorInterface::POSITIVE_FACE:
            f0l = _rheology->getKernelf0Pos(coordsys.get());
            break;
        default:
            PYLITH_COMPONENT_LOGICERROR("Unknown interface face ("<<face<<").");
        } // switch

        kernels.resize(1);
        const EquationPart eqnPart = pylith::feassemble::Integrator::LHS_WEIGHTED;
        kernels[0] = InterfaceResidualKernels("lagrange_multiplier_fault", eqnPart, face, f0l, f1l);
        break;
    } // DYNAMIC_IMEX
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_RETURN(kernels);
} // getInterfaceKernelsResidual


// ------------------------------------------------------------------------------------------------
// Get Jacobian kernels for an interior interface bounding material.
std::vector<pylith::materials::Material::InterfaceJacobianKernels>
pylith::materials::Elasticity::getInterfaceKernelsJacobian(const pylith::topology::Field& solution,
                                                           const pylith::feassemble::IntegratorInterface::FaceEnum face) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getInterfaceKernelsJacobian(solution="<<solution.getName()<<", face="<<face<<")");

    std::vector<InterfaceJacobianKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC:
    case DYNAMIC:
        break;
    case DYNAMIC_IMEX: {
        PetscBdPointJac Jf0ll = nullptr;
        PetscBdPointJac Jf1ll = nullptr;
        PetscBdPointJac Jf2ll = nullptr;
        PetscBdPointJac Jf3ll = nullptr;

        switch (face) {
        case pylith::feassemble::IntegratorInterface::NEGATIVE_FACE:
            Jf0ll = pylith::fekernels::FaultCohesiveKin::Jf0ll_neg;
            break;
        case pylith::feassemble::IntegratorInterface::POSITIVE_FACE:
            Jf0ll = pylith::fekernels::FaultCohesiveKin::Jf0ll_pos;
            break;
        default:
            PYLITH_COMPONENT_LOGICERROR("Unknown interface face ("<<face<<").");
        } // switch

        kernels.resize(1);
        EquationPart eqnPart = pylith::feassemble::Integrator::LHS;
        kernels[0] = InterfaceJacobianKernels("lagrange_multiplier_fault", "lagrange_multiplier_fault", eqnPart, face,
                                              Jf0ll, Jf1ll, Jf2ll, Jf3ll);
        break;
    } // DYNAMIC_IMEX
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    PYLITH_METHOD_RETURN(kernels);
} // getInterfaceKernelsJacobian


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::Elasticity::_updateKernelConstants(const pylith::real dt) {
    assert(_rheology);
    _rheology->updateKernelConstants(&_kernelConstants, dt);
} // _updateKernelConstants


// ------------------------------------------------------------------------------------------------
// Set kernels for residual.
void
pylith::materials::Elasticity::_setKernelsResidual(pylith::feassemble::IntegratorDomain* integrator,
                                                   const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsResidual(integrator="<<integrator<<", solution="<<solution.getName()<<")");

    const std::shared_ptr<spatialdata::geocoords::CoordSys>& coordsys = solution.getMesh().getCoordSys();

    const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
    const int bitGravity = _gravityField ? 0x2 : 0x0;
    const int bitUse = bitBodyForce | bitGravity;

    PetscPointFunc r0 = nullptr;
    switch (bitUse) {
    case 0x1:
        r0 = pylith::fekernels::Elasticity::g0v_bodyforce;
        break;
    case 0x2:
        r0 = pylith::fekernels::Elasticity::g0v_grav;
        break;
    case 0x3:
        r0 = pylith::fekernels::Elasticity::g0v_gravbodyforce;
        break;
    case 0x0:
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for residual kernels.");
    } // switch
    const PetscPointFunc r1 = _rheology->getKernelf1v(coordsys.get());

    std::vector<ResidualKernels> kernels;
    switch (_formulation) {
    case QUASISTATIC: {
        const PetscPointFunc f0u = r0;
        const PetscPointFunc f1u = r1;

        kernels.resize(1);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
        break;
    } // QUASISTATIC
    case DYNAMIC: {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = nullptr;
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = nullptr;

        // Velocity
        const PetscPointFunc f0v = pylith::fekernels::Elasticity::f0v;
        const PetscPointFunc f1v = nullptr;
        const PetscPointFunc g0v = r0;
        const PetscPointFunc g1v = r1;

        kernels.resize(4);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
        kernels[1] = ResidualKernels("displacement", pylith::feassemble::Integrator::RHS, g0u, g1u);
        kernels[2] = ResidualKernels("velocity", pylith::feassemble::Integrator::LHS, f0v, f1v);
        kernels[3] = ResidualKernels("velocity", pylith::feassemble::Integrator::RHS, g0v, g1v);
        break;
    } // DYNAMIC
    case DYNAMIC_IMEX: {
        // Displacement
        const PetscPointFunc f0u = pylith::fekernels::DispVel::f0u;
        const PetscPointFunc f1u = nullptr;
        const PetscPointFunc g0u = pylith::fekernels::DispVel::g0u;
        const PetscPointFunc g1u = nullptr;

        // Velocity
        const PetscPointFunc f0v = pylith::fekernels::DispVel::f0v;
        const PetscPointFunc f1v = nullptr;
        const PetscPointFunc g0v = r0;
        const PetscPointFunc g1v = r1;

        kernels.resize(4);
        kernels[0] = ResidualKernels("displacement", pylith::feassemble::Integrator::LHS, f0u, f1u);
        kernels[1] = ResidualKernels("displacement", pylith::feassemble::Integrator::RHS, g0u, g1u);
        kernels[2] = ResidualKernels("velocity", pylith::feassemble::Integrator::LHS, f0v, f1v);
        kernels[3] = ResidualKernels("velocity", pylith::feassemble::Integrator::RHS, g0v, g1v);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    // Add any MMS body force kernels.
    kernels.insert(kernels.end(), _mmsBodyForceKernels.begin(), _mmsBodyForceKernels.end());

    assert(integrator);
    integrator->setKernelsResidual(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsResidual


// ------------------------------------------------------------------------------------------------
// Set kernels for Jacobian.
void
pylith::materials::Elasticity::_setKernelsJacobian(pylith::feassemble::IntegratorDomain* integrator,
                                                   const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsJacobian(integrator="<<integrator<<", solution="<<solution.getName()<<")");
    assert(integrator);

    const std::shared_ptr<spatialdata::geocoords::CoordSys>& coordsys = solution.getMesh().getCoordSys();

    std::vector<JacobianKernels> kernels;

    switch (_formulation) {
    case QUASISTATIC: {
        const PetscPointJac Jf0uu = nullptr;
        const PetscPointJac Jf1uu = nullptr;
        const PetscPointJac Jf2uu = nullptr;
        const PetscPointJac Jf3uu = _rheology->getKernelJf3vu(coordsys.get());

        integrator->setLHSJacobianTriggers(_rheology->getLHSJacobianTriggers());

        kernels.resize(1);
        const EquationPart equationPart = pylith::feassemble::Integrator::LHS;
        kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        break;
    } // QUASISTATIC
    case DYNAMIC_IMEX: {
        typedef pylith::feassemble::JacobianValues::JacobianKernel ValueKernel;
        std::vector<ValueKernel> valueKernelsJacobian(2);
        std::vector<ValueKernel> valueKernelsPrecond;
        valueKernelsJacobian[0] = ValueKernel("displacement", "displacement", pylith::feassemble::JacobianValues::blockDiag_tshift);
        valueKernelsJacobian[1] = ValueKernel("velocity", "velocity", pylith::feassemble::JacobianValues::blockDiag_tshift);
        integrator->setKernelsJacobian(valueKernelsJacobian, valueKernelsPrecond);

    } // DYNAMIC_IMEX continue with DYNAMIC
    case DYNAMIC: {
        const PetscPointJac Jf0uu = pylith::fekernels::DispVel::Jf0uu_stshift;
        const PetscPointJac Jf1uu = nullptr;
        const PetscPointJac Jf2uu = nullptr;
        const PetscPointJac Jf3uu = nullptr;

        const PetscPointJac Jf0vv = pylith::fekernels::Elasticity::Jf0vv;
        const PetscPointJac Jf1vv = nullptr;
        const PetscPointJac Jf2vv = nullptr;
        const PetscPointJac Jf3vv = nullptr;

        integrator->setLHSJacobianTriggers(pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE);

        kernels.resize(2);
        const EquationPart equationPart = pylith::feassemble::Integrator::LHS_LUMPED_INV;
        kernels[0] = JacobianKernels("displacement", "displacement", equationPart, Jf0uu, Jf1uu, Jf2uu, Jf3uu);
        kernels[1] = JacobianKernels("velocity", "velocity", equationPart, Jf0vv, Jf1vv, Jf2vv, Jf3vv);
        break;
    } // DYNAMIC
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown formulation for equations (" << _formulation << ").");
    } // switch

    integrator->setKernelsJacobian(kernels, solution);

    PYLITH_METHOD_END;
} // _setKernelsJacobian


// ------------------------------------------------------------------------------------------------
// Set kernels for computing updated state variables in auxiliary field.
void
pylith::materials::Elasticity::_setKernelsUpdateStateVars(pylith::feassemble::IntegratorDomain* integrator,
                                                          const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsUpdateStateVars(integrator="<<integrator<<", solution="<<solution.getName()<<")");

    const std::shared_ptr<spatialdata::geocoords::CoordSys>& coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels;
    _rheology->addKernelsUpdateStateVars(&kernels, coordsys.get());

    integrator->setKernelsUpdateStateVars(kernels);

    PYLITH_METHOD_END;
} // _setKernelsUpdateStateVars


// ------------------------------------------------------------------------------------------------
// Set kernels for computing derived field.
void
pylith::materials::Elasticity::_setKernelsDerivedField(pylith::feassemble::IntegratorDomain* integrator,
                                                       const topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_setKernelsDerivedField(integrator="<<integrator<<", solution="<<solution.getName()<<")");

    const std::shared_ptr<spatialdata::geocoords::CoordSys>& coordsys = solution.getMesh().getCoordSys();
    assert(coordsys);

    std::vector<ProjectKernels> kernels(2);
    kernels[0] = ProjectKernels("cauchy_stress", _rheology->getKernelCauchyStressVector(coordsys.get()));

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc strainKernel =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
        nullptr;
    kernels[1] = ProjectKernels("cauchy_strain", strainKernel);

    assert(integrator);
    integrator->setKernelsDerivedField(kernels);

    PYLITH_METHOD_END;
} // _setKernelsDerivedField


// End of file
