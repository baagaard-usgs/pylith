// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/materials/IsotropicLinearElasticity.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryElastic.hh" // USES AuxiliaryFactoryElastic
#include "pylith/fekernels/IsotropicLinearElasticity.hh" // USES IsotropicLinearElasticity kernels
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <typeinfo> // USES typeid()

namespace pylith {
    namespace materials {
        class _IsotropicLinearElasticity {
public:

            /** Validate auxiliary field.
             *
             * @param[in] auxiliaryField Auxiliary field to validate.
             */
            static void validateAuxiliaryField(const pylith::topology::Field* auxiliaryField);

        };
    }
}

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
// Get validator for auxiliary field.
pylith::feassemble::AuxiliaryFactory::validatorfn_type
pylith::materials::IsotropicLinearElasticity::getAuxiliaryValidator(void) {
    return _IsotropicLinearElasticity::validateAuxiliaryField;
}


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


// ------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelResidualStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1v =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1v_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f1v_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f1v_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(f1v);
} // getKernelf1v


// ------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearElasticity::getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJacobianElasticConstants(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3vu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::Jf3vu_infinitesimalStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::Jf3vu_infinitesimalStrain :
        NULL;

    PYLITH_METHOD_RETURN(Jf3vu);
} // getKernelJacobianElasticConstants


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelf0Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Neg(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_neg_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_neg_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_neg_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_neg_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelf0Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_pos_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_pos_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::f0l_pos_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::f0l_pos_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearElasticity::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelDerivedCauchyStress(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_infinitesimalStrain_asVector :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticity3D::cauchyStress_infinitesimalStrain_refState_asVector :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearElasticityPlaneStrain::cauchyStress_infinitesimalStrain_refState_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelDerivedCauchyStress


// ------------------------------------------------------------------------------------------------
// Validate auxiliary field.
void
pylith::materials::_IsotropicLinearElasticity::validateAuxiliaryField(const pylith::topology::Field* auxiliaryField) {
    assert(auxiliaryField);
    const PylithReal tolerance = pylith::feassemble::AuxiliaryFactory::SCALE_TOLERANCE;

    const PylithReal minValue = 1.0 / tolerance;
    const PylithReal maxValue = tolerance;

    const PetscInt i_shearModulus = auxiliaryField->getSubfieldInfo("shear_modulus").index;
    const PetscInt i_bulkModulus = auxiliaryField->getSubfieldInfo("bulk_modulus").index;
    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField);
    const PetscScalar* auxiliaryArray = auxiliaryVisitor.localArray();

    PetscInt pStart = 0, pEnd = 0;
    PetscErrorCode err = PETSC_SUCCESS;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    PetscReal minShearModulus = std::numeric_limits<float>::max();
    PetscReal maxShearModulus = std::numeric_limits<float>::min();
    PetscReal minBulkModulus = std::numeric_limits<float>::max();
    PetscReal maxBulkModulus = std::numeric_limits<float>::min();
    for (PetscInt point = pStart; point < pEnd; ++point) {
        PetscInt dof = auxiliaryVisitor.sectionSubfieldDof(i_shearModulus, point);
        if (dof) {
            const PetscInt offset = auxiliaryVisitor.sectionSubfieldOffset(i_shearModulus, point);
            const PetscScalar shearModulus = auxiliaryArray[offset];
            minShearModulus = std::min(minShearModulus, shearModulus);
            maxShearModulus = std::max(maxShearModulus, shearModulus);
        } // if

        dof = auxiliaryVisitor.sectionSubfieldDof(i_bulkModulus, point);
        if (dof) {
            const PetscInt offset = auxiliaryVisitor.sectionSubfieldOffset(i_bulkModulus, point);
            const PetscScalar bulkModulus = auxiliaryArray[offset];
            minBulkModulus = std::min(minBulkModulus, bulkModulus);
            maxBulkModulus = std::max(maxBulkModulus, bulkModulus);
        } // if
    } // for

    const bool shearModulusFailures = minShearModulus < minValue || maxShearModulus > maxValue;
    const bool bulkModulusFailures = minBulkModulus < minValue || maxBulkModulus > maxValue;
    if (shearModulusFailures || bulkModulusFailures) {
        std::ostringstream msg;
        msg << "Auxiliary field '" << auxiliaryField->getLabel() << "' failed nondimensioanlization check.\n";
        if (shearModulusFailures) {
            msg << "Found nondimensional shear modulus value(s) outside of range ["<<minValue<<", "<<maxValue<<"]. "
                << "Minimum value: " << minShearModulus << ", maximum value: " << maxShearModulus << "\n";
        } // if
        if (bulkModulusFailures) {
            msg << "Found nondimensional bulk modulus value(s) outside of range ["<<minValue<<", "<<maxValue<<"]. "
                << "Minimum value: " << minBulkModulus << ", maximum value: " << maxBulkModulus << "\n";
        } // if
        throw std::runtime_error(msg.str());
    } // if
}


// End of file
