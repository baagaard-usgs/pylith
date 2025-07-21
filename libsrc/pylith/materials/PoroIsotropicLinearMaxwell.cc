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

#include "pylith/materials/PoroIsotropicLinearMaxwell.hh" // implementation of object methods
#include "pylith/materials/AuxiliaryFactoryPoroelastic.hh" // USES AuxiliaryFactory

// #include "pylith/fekernels/IsotropicLinearPoroelasticity.hh" // USES IsotropicLinearIncompElasticity kernels
#include "pylith/fekernels/PoroIsotropicLinearMaxwell.hh" // USES IsotropicLinearIncompElasticity kernels
#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()
#include "PoroIsotropicLinearMaxwell.hh"

// ---------------------------------------------------------------------------------------------------------------------
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::PoroIsotropicLinearMaxwell::PoroIsotropicLinearMaxwell(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryPoroelastic),
    _useReferenceState(false),
    _useTensorPermeability(false) {
    pylith::utils::PyreComponent::setName("poroisotropiclinearmaxwell");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::PoroIsotropicLinearMaxwell::~PoroIsotropicLinearMaxwell(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::PoroIsotropicLinearMaxwell::deallocate(void) {
    RheologyPoroelasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
void
pylith::materials::PoroIsotropicLinearMaxwell::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and strain?
bool
pylith::materials::PoroIsotropicLinearMaxwell::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ---------------------------------------------------------------------------------------------------------------------
// Use full tensor permeability?
void
pylith::materials::PoroIsotropicLinearMaxwell::useTensorPermeability(const bool value) {
    PYLITH_COMPONENT_DEBUG("useTensorPermeability="<<value<<")");

    _useTensorPermeability = value;
} // useTensorPermeability


// ---------------------------------------------------------------------------------------------------------------------
// Use full tensor permeability?
bool
pylith::materials::PoroIsotropicLinearMaxwell::useTensorPermeability(void) const {
    return _useTensorPermeability;
} // useTensorPermeability


// ---------------------------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryPoroelastic*
pylith::materials::PoroIsotropicLinearMaxwell::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ---------------------------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::PoroIsotropicLinearMaxwell::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress(); // numA - 10
        _auxiliaryFactory->addReferenceStrain(); // numA - 9
    } // if
    _auxiliaryFactory->addShearModulus(); // Shear Modulus, G, numA - 8
    _auxiliaryFactory->addDrainedBulkModulus(); // K_d, numA - 7
    _auxiliaryFactory->addBiotCoefficient(); // alpha, numA - 6
    _auxiliaryFactory->addBiotModulus(); // M, numA - 5
    _auxiliaryFactory->addMaxwellTime(); // numA - 4
    _auxiliaryFactory->addViscousStrain(); // numA - 3
    _auxiliaryFactory->addTotalStrain(); // numA - 2
    if (_useTensorPermeability) {
        _auxiliaryFactory->addTensorPermeability(); // k, numA - 1
    } else {
        _auxiliaryFactory->addIsotropicPermeability(); // k, numA - 1
    }

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// =============================== LHS =========================================

// ---------------------------------------------------------------------------------------------------------------------
// Select implicit f0p function.
PetscPointFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelf0p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                     const bool _useBodyForce,
                                                                     const bool _gravityField,
                                                                     const bool _useSourceDensity) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0p="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitBodyForce = _useBodyForce ? 0x1 : 0x0;
    const int bitGravity = _gravityField ? 0x2 : 0x0;
    const int bitSourceDensity = _useSourceDensity ? 0x4 : 0x0;
    const int bitUse = bitBodyForce | bitGravity | bitSourceDensity;

    PetscPointFn* f0p = NULL;
    switch (bitUse) {
    case 0x0:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x1:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x2:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x4:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit_source :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit_source :
              NULL; // aOff for sourceDensity is 3
        break;
    case 0x3:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit :
              NULL;
        break;
    case 0x5:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit_source_body :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit_source_body :
              NULL; // aOff for sourceDensity is 4
        break;
    case 0x6:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit_source_grav :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit_source_grav :
              NULL; // aOff for sourceDensity is 4
        break;
    case 0x7:
        f0p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f0p_implicit_source_grav_body :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f0p_implicit_source_grav_body :
              NULL; // aOff for sourceDensity is 5
        break;
    default:
        PYLITH_COMPONENT_LOGICERROR("Unknown case (bitUse=" << bitUse << ") for Poroelasticity LHS residual kernels.");
    } // switch
    PYLITH_METHOD_RETURN(f0p);
} // getKernelf0p_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelf1u_implicit(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1u(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitReferenceState = _useReferenceState ? 0x1 : 0x0;
    const int bitUse = bitReferenceState;

    PetscPointFn* f1u = NULL;
    switch (bitUse) {
    case 0x0:
        f1u = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1u :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1u :
              NULL;
        break;
    case 0x1:
        f1u = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1u_refstate :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1u_refstate :
              NULL;
        break;
    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for getKernelf1u_implicit.");
        throw std::logic_error("Unknown combination of flags.");
    } // switch
    PYLITH_METHOD_RETURN(f1u);
} // getKernelf1u_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get darcy velocity kernel
PetscPointFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelf1p_implicit(const spatialdata::geocoords::CoordSys* coordsys,
                                                                     const bool _useBodyForce,
                                                                     const bool _gravityField) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1p_implicit="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const int bitTensorPermeability = _useTensorPermeability ? 0x1 : 0x0;
    const int bitBodyForce = _useBodyForce ? 0x2 : 0x0;
    const int bitGravityField = _gravityField ? 0x4 : 0x0;
    const int bitUse = bitTensorPermeability | bitBodyForce | bitGravityField;

    PetscPointFn* f1p = NULL;
    switch (bitUse) {
    case 0x0:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p :
              NULL;
        break;
    case 0x1:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p_tensor_permeability :
              NULL;
        break;
    case 0x2:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p_body :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p_body :
              NULL;
        break;
    case 0x3:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p_body_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p_body_tensor_permeability :
              NULL;
        break;
    case 0x4:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p_gravity :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p_gravity :
              NULL;
        break;
    case 0x5:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p_gravity_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p_gravity_tensor_permeability :
              NULL;
        break;
    case 0x6:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p_body_gravity :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p_body_gravity :
              NULL;
        break;
    case 0x7:
        f1p = (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::f1p_body_gravity_tensor_permeability :
              (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::f1p_body_gravity_tensor_permeability :
              NULL;
        break;

    default:
        PYLITH_COMPONENT_ERROR("Unknown combination of flags for  _useTensorPermeability="<<_useTensorPermeability<<", _useBodyForce="<<_useBodyForce<<", _gravityField="<<_gravityField<<").");
        throw std::logic_error("Unknown combination of flags.");
    } // switch
    PYLITH_METHOD_RETURN(f1p);
} // getKernelf1p_implicit


// ---------------------------------------------------------------------------------------------------------------------
// Get poroelastic constants kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelJf3uu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3uu(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::Jf3uu :
        (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::Jf3uu :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
    // PYLITH_METHOD_RETURN(NULL);

} // getKernelLHSJacobianElasticConstants


// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelJf2up(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2up(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf2up =
        (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::Jf2up :
        (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::Jf2up :
        NULL;

    PYLITH_METHOD_RETURN(Jf2up);
    // PYLITH_METHOD_RETURN(NULL);

} // getKernelJf2up


// ---------------------------------------------------------------------------------------------------------------------
// Get lambda kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelJf2ue(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf2ue(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf2ue =
        (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::Jf2ue :
        (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::Jf2ue :
        NULL;

    PYLITH_METHOD_RETURN(Jf2ue);
    // PYLITH_METHOD_RETURN(NULL);
} // getKernelJf2ue


// ---------------------------------------------------------------------------------------------------------------------
// Get Specific storage kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJacFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelJf0pp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pp(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf0pp =
        (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::Jf0pp :
        (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::Jf0pp :
        NULL;

    PYLITH_METHOD_RETURN(Jf0pp);
    // PYLITH_METHOD_RETURN(NULL);
} // getKernelJf0pp


// ---------------------------------------------------------------------------------------------------------------------
// Get Darcy Conductivity kernel for LHS Jacobian
PetscPointJacFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelJf3pp(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3pp(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf3pp =
        (!_useTensorPermeability && 3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::Jf3pp :
        (!_useTensorPermeability && 2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::Jf3pp :
        (_useTensorPermeability && 3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::Jf3pp_tensor_permeability :
        (_useTensorPermeability && 2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::Jf3pp_tensor_permeability :
        NULL;

    PYLITH_METHOD_RETURN(Jf3pp);
    // PYLITH_METHOD_RETURN(NULL);

} // getKerneJf3pp


// ---------------------------------------------------------------------------------------------------------------------
// Get biot coefficient kernel for LHS Jacobian F(t,s, \dot{s}).
PetscPointJacFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelJf0pe(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf0pe(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJacFn* Jf0pe =
        (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::Jf0pe :
        (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::Jf0pe :
        NULL;

    PYLITH_METHOD_RETURN(Jf0pe);
    // PYLITH_METHOD_RETURN(NULL);
} // getKernelJf0pe


// =========================== DERIVED FIELDS ==================================

// ---------------------------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelCauchyStressVector(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* kernel =
        (!_useReferenceState && 3 == spaceDim) ?  pylith::fekernels::PoroIsotropicLinearMaxwell3D::cauchyStress_infinitesimalStrain_asVector :
        (!_useReferenceState && 2 == spaceDim) ?  pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        (_useReferenceState && 3 == spaceDim) ?  pylith::fekernels::PoroIsotropicLinearMaxwell3D::cauchyStress_infinitesimalStrain_refState_asVector :
        (_useReferenceState && 2 == spaceDim) ?  pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::cauchyStress_infinitesimalStrain_refState_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector


// ---------------------------------------------------------------------------------------------------------------------
// Get water content kernel for derived field.
PetscPointFn*
pylith::materials::PoroIsotropicLinearMaxwell::getKernelWaterContent(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelWaterContent(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* kernel =
        (3 == spaceDim) ?  pylith::fekernels::PoroIsotropicLinearMaxwell3D::waterContent_asScalar :
        (2 == spaceDim) ?  pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::waterContent_asScalar :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelWaterContent


// ---------------------------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::PoroIsotropicLinearMaxwell::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                     const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<")");

    assert(kernelConstants);

    if (1 != kernelConstants->size()) { kernelConstants->resize(1);}
    (*kernelConstants)[0] = dt;

    PYLITH_METHOD_END;
} // updateKernelConstants


// ---------------------------------------------------------------------------------------------------------------------
// Add kernels for updating state variables, implicit.
void
pylith::materials::PoroIsotropicLinearMaxwell::addKernelsUpdateStateVarsImplicit(std::vector<ProjectKernels>* kernels,
                                                                                 const spatialdata::geocoords::CoordSys* coordsys,
                                                                                 const bool _useStateVars) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVarsImplicit(kernels="<<kernels<<", coordsys="<<coordsys<<")");
    if (_useStateVars) {
        const int spaceDim = coordsys->getSpaceDim();

        const PetscPointFn* funcPorosity =
            (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::updatePorosityImplicit :
            (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::updatePorosityImplicit :
            NULL;
        const PetscPointFn* funcViscousStrain =
            (3 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwell3D::viscousStrain_infinitesimalStrain_asVector :
            (2 == spaceDim) ? pylith::fekernels::PoroIsotropicLinearMaxwellPlaneStrain::viscousStrain_infinitesimalStrain_asVector :
            NULL;
        const PetscPointFn* funcTotalStrain =
            (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
            (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
            NULL;

        assert(kernels);
        size_t prevNumKernels = kernels->size();
        kernels->resize(prevNumKernels + 1);
        (*kernels)[prevNumKernels+0] = ProjectKernels("porosity", funcPorosity);
        (*kernels)[prevNumKernels+1] = ProjectKernels("viscous_strain", funcViscousStrain);
        (*kernels)[prevNumKernels+2] = ProjectKernels("total_strain", funcTotalStrain);
    }

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVarsImplicit


// End of file
