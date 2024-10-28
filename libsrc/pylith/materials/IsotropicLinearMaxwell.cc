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

#include "pylith/materials/IsotropicLinearMaxwell.hh" // implementation of object methods

#include "pylith/materials/AuxiliaryFactoryViscoelastic.hh" // USES AuxiliaryFactoryViscoelastic
#include "pylith/fekernels/IsotropicLinearMaxwell.hh" // USES IsotropicLinearMaxwell kernels
#include "pylith/feassemble/Integrator.hh" // USES Integrator
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <typeinfo> // USES typeid()

namespace pylith {
    namespace materials {
        class _IsotropicLinearMaxwell {
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
typedef pylith::feassemble::IntegratorDomain::ProjectKernels ProjectKernels;

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::IsotropicLinearMaxwell::IsotropicLinearMaxwell(void) :
    _auxiliaryFactory(new pylith::materials::AuxiliaryFactoryViscoelastic),
    _useReferenceState(false) {
    _lhsJacobianTriggers = pylith::feassemble::Integrator::NEW_JACOBIAN_TIME_STEP_CHANGE;
    pylith::utils::PyreComponent::setName("isotropiclinearmaxwell");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::IsotropicLinearMaxwell::~IsotropicLinearMaxwell(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::IsotropicLinearMaxwell::deallocate(void) {
    RheologyElasticity::deallocate();

    delete _auxiliaryFactory;_auxiliaryFactory = NULL;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
void
pylith::materials::IsotropicLinearMaxwell::useReferenceState(const bool value) {
    PYLITH_COMPONENT_DEBUG("useReferenceState="<<value<<")");

    _useReferenceState = value;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Use reference stress and strain in computation of stress and
// strain?
bool
pylith::materials::IsotropicLinearMaxwell::useReferenceState(void) const {
    return _useReferenceState;
} // useReferenceState


// ------------------------------------------------------------------------------------------------
// Get auxiliary factory associated with physics.
pylith::materials::AuxiliaryFactoryElasticity*
pylith::materials::IsotropicLinearMaxwell::getAuxiliaryFactory(void) {
    return _auxiliaryFactory;
} // getAuxiliaryFactory


// ------------------------------------------------------------------------------------------------
// Get validator for auxiliary field.
pylith::feassemble::AuxiliaryFactory::validatorfn_type
pylith::materials::IsotropicLinearMaxwell::getAuxiliaryValidator(void) {
    return _IsotropicLinearMaxwell::validateAuxiliaryField;
}


// ------------------------------------------------------------------------------------------------
// Add rheology subfields to auxiliary field.
void
pylith::materials::IsotropicLinearMaxwell::addAuxiliarySubfields(void) {
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
    _auxiliaryFactory->addMaxwellTime();
    _auxiliaryFactory->addViscousStrain();
    _auxiliaryFactory->addTotalStrain();

    PYLITH_METHOD_END;
} // addAuxiliarySubfields


// ------------------------------------------------------------------------------------------------
// Get stress kernel for LHS residual, F(t,s,\dot{s}).
PetscPointFunc
pylith::materials::IsotropicLinearMaxwell::getKernelf1v(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf1v(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc f1u =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f1v_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f1v_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f1v_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(f1u);
} // getKernelf1v


// ------------------------------------------------------------------------------------------------
// Get elastic constants kernel for LHS Jacobian F(t,s,\dot{s}).
PetscPointJac
pylith::materials::IsotropicLinearMaxwell::getKernelJf3vu(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelJf3vu(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointJac Jf3uu =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::Jf3vu_infinitesimalStrain :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::Jf3vu_infinitesimalStrain :
        NULL;

    PYLITH_METHOD_RETURN(Jf3uu);
} // getKernelJf3vu


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for negative fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearMaxwell::getKernelf0Neg(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Neg(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_neg_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_neg_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_neg_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_neg_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get f0 kernel for LHS interface residual, F(t,s), for positive fault face.
PetscBdPointFunc
pylith::materials::IsotropicLinearMaxwell::getKernelf0Pos(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelf0Pos(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscBdPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_pos_infinitesimalStrain :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_pos_infinitesimalStrain :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::f0l_pos_infinitesimalStrain_refState :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::f0l_pos_infinitesimalStrain_refState :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
}


// ------------------------------------------------------------------------------------------------
// Get stress kernel for derived field.
PetscPointFunc
pylith::materials::IsotropicLinearMaxwell::getKernelCauchyStressVector(const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("getKernelCauchyStressVector(coordsys="<<typeid(coordsys).name()<<")");

    const int spaceDim = coordsys->getSpaceDim();
    PetscPointFunc kernel =
        (!_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress_infinitesimalStrain_asVector :
        (!_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress_infinitesimalStrain_asVector :
        (_useReferenceState && 3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::cauchyStress_infinitesimalStrain_refState_asVector :
        (_useReferenceState && 2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::cauchyStress_infinitesimalStrain_refState_asVector :
        NULL;

    PYLITH_METHOD_RETURN(kernel);
} // getKernelCauchyStressVector


// ------------------------------------------------------------------------------------------------
// Update kernel constants.
void
pylith::materials::IsotropicLinearMaxwell::updateKernelConstants(pylith::real_array* kernelConstants,
                                                                 const PylithReal dt) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("updateKernelConstants(kernelConstants"<<kernelConstants<<", dt="<<dt<<")");

    assert(kernelConstants);

    if (1 != kernelConstants->size()) { kernelConstants->resize(1);}
    (*kernelConstants)[0] = dt;

    PYLITH_METHOD_END;
} // updateKernelConstants


// ------------------------------------------------------------------------------------------------
// Add kernels for updating state variables.
void
pylith::materials::IsotropicLinearMaxwell::addKernelsUpdateStateVars(std::vector<ProjectKernels>* kernels,
                                                                     const spatialdata::geocoords::CoordSys* coordsys) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("addKernelsUpdateStateVars(kernels="<<kernels<<", coordsys="<<coordsys<<")");

    const int spaceDim = coordsys->getSpaceDim();
    const PetscPointFunc funcViscousStrain =
        (3 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwell3D::viscousStrain_infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::IsotropicLinearMaxwellPlaneStrain::viscousStrain_infinitesimalStrain_asVector :
        NULL;
    const PetscPointFunc funcTotalStrain =
        (3 == spaceDim) ? pylith::fekernels::Elasticity3D::infinitesimalStrain_asVector :
        (2 == spaceDim) ? pylith::fekernels::ElasticityPlaneStrain::infinitesimalStrain_asVector :
        NULL;

    assert(kernels);
    size_t prevNumKernels = kernels->size();
    kernels->resize(prevNumKernels + 2);
    (*kernels)[prevNumKernels+0] = ProjectKernels("viscous_strain", funcViscousStrain);
    (*kernels)[prevNumKernels+1] = ProjectKernels("total_strain", funcTotalStrain);

    PYLITH_METHOD_END;
} // addKernelsUpdateStateVars


// ------------------------------------------------------------------------------------------------
// Validate auxiliary field.
void
pylith::materials::_IsotropicLinearMaxwell::validateAuxiliaryField(const pylith::topology::Field* auxiliaryField) {
    assert(auxiliaryField);
    const PylithReal tolerance = pylith::feassemble::AuxiliaryFactory::SCALE_TOLERANCE;

    const PylithReal minValue = 1.0 / tolerance;
    const PylithReal maxValue = tolerance;

    const PetscInt i_shearModulus = auxiliaryField->getSubfieldInfo("shear_modulus").index;
    const PetscInt i_bulkModulus = auxiliaryField->getSubfieldInfo("bulk_modulus").index;
    const PetscInt i_maxwellTime = auxiliaryField->getSubfieldInfo("maxwell_time").index;
    pylith::topology::VecVisitorMesh auxiliaryVisitor(*auxiliaryField);
    const PetscScalar* auxiliaryArray = auxiliaryVisitor.localArray();

    PetscInt pStart = 0, pEnd = 0;
    PetscErrorCode err = PETSC_SUCCESS;
    err = PetscSectionGetChart(auxiliaryField->getLocalSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

    PetscReal minShearModulus = std::numeric_limits<float>::max();
    PetscReal maxShearModulus = std::numeric_limits<float>::min();
    PetscReal minBulkModulus = std::numeric_limits<float>::max();
    PetscReal maxBulkModulus = std::numeric_limits<float>::min();
    PetscReal minMaxwellTime = std::numeric_limits<float>::max();
    PetscReal maxMaxwellTime = std::numeric_limits<float>::min();
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

        dof = auxiliaryVisitor.sectionSubfieldDof(i_maxwellTime, point);
        if (dof) {
            const PetscInt offset = auxiliaryVisitor.sectionSubfieldOffset(i_maxwellTime, point);
            const PetscScalar maxwellTime = auxiliaryArray[offset];
            minMaxwellTime = std::min(minMaxwellTime, maxwellTime);
            maxMaxwellTime = std::max(maxMaxwellTime, maxwellTime);
        } // if
    } // for

    const bool shearModulusFailures = minShearModulus < minValue || maxShearModulus > maxValue;
    const bool bulkModulusFailures = minBulkModulus < minValue || maxBulkModulus > maxValue;
    const bool maxwellTimeFailures = minMaxwellTime < minValue || maxMaxwellTime > maxValue;
    if (shearModulusFailures || bulkModulusFailures || maxwellTimeFailures) {
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
        if (maxwellTimeFailures) {
            msg << "Found nondimensional Maxwell time value(s) outside of range ["<<minValue<<", "<<maxValue<<"]. "
                << "Minimum value: " << minMaxwellTime << ", maximum value: " << maxMaxwellTime << "\n";
        } // if
        throw std::runtime_error(msg.str());
    } // if
}


// End of file
