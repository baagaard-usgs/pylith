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

#include "pylith/materials/IsotropicLinearElasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearElasticity::IsotropicLinearElasticity(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryElastic),
    _useReferenceState(false) {
    pylith::utils::PyreComponent::setName("isotropiclinearelasticity");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearElasticity::~IsotropicLinearElasticity(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearElasticity::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearElasticity::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearElasticity::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearElasticity::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearElasticity::addAuxiliarySubfields(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addAuxiliarySubfields(void)");

    // :ATTENTION: The order for adding subfields must match the order of the auxiliary fields in the point-wise
    // functions (kernels).

    if (_useReferenceState) {
        _auxiliaryFactory->addReferenceStress();
        _auxiliaryFactory->addReferenceStrain();
    } // if
    _auxiliaryFactory->addShearModulus();
    _auxiliaryFactory->addBulkModulus();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


#include "pylith/fekernels//pde/elasticity/SolutionLayout.hh"
#include "pylith/fekernels//pde/elasticity/IsotropicLinearLayout.hh"
#include "pylith/fekernels/momentum/pde/DivStress.hh"
#include "pylith/fekernels/momentum/strain/InfinitesimalStrain.hh"
#include "pylith/fekernels/momentum/stress/elasticity/IsotropicLinear.hh"
// ------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFn*
pylith::materials::IsotropicLinearElasticity::getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    // :TODO: Change name to getResidualKernelMomentumDivStress

    // :TODO: Make these arguments to function.
    constexpr pylith::fekernels::pde::elasticity::SolutionFlags solnFlags = pylith::fekernels::pde::elasticity::DEFAULT;
    constexpr pylith::fekernels::momentum::MomentumFlags momentumFlags = pylith::fekernels::momentum::DEFAULT;

    const int spaceDim = coordsys->getSpaceDim();
    using SolutionLayout = pylith::fekernels::pde::elasticity::SolutionLayout<solnFlags>;
    PetscPointFn* f1v = nullptr;
    if (_useReferenceState) {
        constexpr pylith::fekernels::pde::elasticity::IsotropicLinearFlags auxFlags = pylith::fekernels::pde::elasticity::ISOTROPIC_LINEAR_DEFAULT | pylith::fekernels::pde::elasticity::ISOTROPIC_LINEAR_REFERENCE_STRAIN | pylith::fekernels::pde::elasticity::ISOTROPIC_LINEAR_REFERENCE_STRESS;
        using AuxiliaryLayout = pylith::fekernels::pde::elasticity::IsotropicLinearLayout<momentumFlags, auxFlags>;

        if (2 == spaceDim) {
            constexpr size_t dim = 2;
            using StrainModel = pylith::fekernels::momentum::InfinitesimalStrain<dim, SolutionLayout>;
            using StressModel = pylith::fekernels::momentum::stress::elasticity::IsotropicLinear<dim, AuxiliaryLayout>;
            f1v = pylith::fekernels::momentum::DivStress<dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout>::f1;
        } else if (3 == spaceDim) {
            constexpr size_t dim = 3;
            using StrainModel = pylith::fekernels::momentum::InfinitesimalStrain<dim, SolutionLayout>;
            using StressModel = pylith::fekernels::momentum::stress::elasticity::IsotropicLinear<dim, AuxiliaryLayout>;
            f1v = pylith::fekernels::momentum::DivStress<dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout>::f1;
        } // if/else

    } else {
        constexpr pylith::fekernels::pde::elasticity::IsotropicLinearFlags auxFlags = pylith::fekernels::pde::elasticity::ISOTROPIC_LINEAR_DEFAULT;
        using AuxiliaryLayout = pylith::fekernels::pde::elasticity::IsotropicLinearLayout<momentumFlags, auxFlags>;

        if (2 == spaceDim) {
            constexpr size_t dim = 2;
            using StrainModel = pylith::fekernels::momentum::InfinitesimalStrain<dim, SolutionLayout>;
            using StressModel = pylith::fekernels::momentum::stress::elasticity::IsotropicLinear<dim, AuxiliaryLayout>;
            f1v = pylith::fekernels::momentum::DivStress<dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout>::f1;
        } else if (3 == spaceDim) {
            constexpr size_t dim = 3;
            using StrainModel = pylith::fekernels::momentum::InfinitesimalStrain<dim, SolutionLayout>;
            using StressModel = pylith::fekernels::momentum::stress::elasticity::IsotropicLinear<dim, AuxiliaryLayout>;
            f1v = pylith::fekernels::momentum::DivStress<dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout>::f1;
        } // if/else
    } // if/else

    PYLITH_METHOD_RETURN(f1v);
} // getKernelf1v


// ------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJacFn*
pylith::materials::IsotropicLinearElasticity::getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    // :TODO: Change name to getJacobianKernelMomentumDivStress

    // :TODO: Make these arguments to function.
    constexpr pylith::fekernels::pde::elasticity::SolutionFlags solnFlags = pylith::fekernels::pde::elasticity::DEFAULT;
    constexpr pylith::fekernels::momentum::MomentumFlags momentumFlags = pylith::fekernels::momentum::DEFAULT;

    const size_t spaceDim = coordsys->getSpaceDim();
    constexpr pylith::fekernels::pde::elasticity::IsotropicLinearFlags auxFlags = pylith::fekernels::pde::elasticity::ISOTROPIC_LINEAR_DEFAULT;
    using SolutionLayout = pylith::fekernels::pde::elasticity::SolutionLayout<solnFlags>;
    using AuxiliaryLayout = pylith::fekernels::pde::elasticity::IsotropicLinearLayout<momentumFlags, auxFlags>;

    PetscPointJacFn* Jf3vu = nullptr;
    if (2 == spaceDim) {
        constexpr size_t dim = 2;
        using StrainModel = pylith::fekernels::momentum::InfinitesimalStrain<dim, SolutionLayout>;
        using StressModel = pylith::fekernels::momentum::stress::elasticity::IsotropicLinear<dim, AuxiliaryLayout>;
        Jf3vu = pylith::fekernels::momentum::DivStress<dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout>::Jf3uu;
    } else if (3 == spaceDim) {
        constexpr size_t dim = 3;
        using StrainModel = pylith::fekernels::momentum::InfinitesimalStrain<dim, SolutionLayout>;
        using StressModel = pylith::fekernels::momentum::stress::elasticity::IsotropicLinear<dim, AuxiliaryLayout>;
        Jf3vu = pylith::fekernels::momentum::DivStress<dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout>::Jf3uu;
    } // if/else

    PYLITH_METHOD_RETURN(Jf3vu);
} // getKernelJacobianElasticConstants


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFn*
pylith::materials::IsotropicLinearElasticity::getKernelf0Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Neg(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFn* kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_neg_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_neg_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_neg_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_neg_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFn*
pylith::materials::IsotropicLinearElasticity::getKernelf0Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFn* kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_pos_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_pos_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_pos_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_pos_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFn*
pylith::materials::IsotropicLinearElasticity::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFn* kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_infinitesimalStrain_asVector :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_infinitesimalStrain_refState_asVector :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_infinitesimalStrain_refState_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// End of file
