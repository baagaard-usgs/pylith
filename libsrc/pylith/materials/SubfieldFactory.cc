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

#include "pylith/materials/SubfieldFactory.hh" // implementation of object methods

#include "pylith/materials/Query.hh" // USES Query

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // HOLDSA FieldQuery

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace materials {
        class _SubfieldFactory {
public:

            ///< Names of field components for 2D tensor.
            static const char* componentsVector[3];

            ///< Names of field components for tensor.
            static const char* componentsTensor[6];

            static const char* genericComponent; ///< Name of generic PyLith component

            /** Set names of vector field components.
             *
             * @pre The number of components must be set in the description.
             *
             * @param[inout] description Description in which to set names.
             */
            static
            void setVectorComponentNames(pylith::topology::FieldBase::Description* description);

            /** Set names of tensor field components.
             *
             * @pre The number of components must be set in the description.
             *
             * @param[inout] description Description in which to set names.
             */
            static
            void setTensorComponentNames(pylith::topology::FieldBase::Description* description);

            static
            pylith::topology::FieldBase::Description getDensity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getBodyForce(const spatialdata::units::Nondimensional& normalizer,
                                                                  const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getGravityField(const spatialdata::units::Nondimensional& normalizer,
                                                                     const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getShearModulus(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getBulkModulus(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getReferenceStress(const spatialdata::units::Nondimensional& normalizer,
                                                                        const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getReferenceStrain(const spatialdata::units::Nondimensional& normalizer,
                                                                        const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getMaxwellTime(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getGenMaxwellMaxwellTime(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getGenMaxwellShearModulusRatio(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getPowerlawReferenceStrainRate(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getPowerlawReferenceStress(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getPowerlawExponent(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTotalStrain(const spatialdata::units::Nondimensional& normalizer,
                                                                    const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getDeviatoricStress(const spatialdata::units::Nondimensional& normalizer,
                                                                         const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getViscousStrain(const spatialdata::units::Nondimensional& normalizer,
                                                                      const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getGenMaxwellViscousStrain(const spatialdata::units::Nondimensional& normalizer,
                                                                                const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getIsotropicPermeability(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTensorPermeability(const spatialdata::units::Nondimensional& normalizer,
                                                                           const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getDrainedBulkModulus(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getBiotCoefficient(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getBiotModulus(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getPorosity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getSolidDensity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getFluidDensity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getFluidViscosity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getSourceDensity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getCauchyStress(const spatialdata::units::Nondimensional& normalizer,
                                                                     const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getCauchyStrain(const spatialdata::units::Nondimensional& normalizer,
                                                                     const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getBulkDensity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getWateContent(const spatialdata::units::Nondimensional& normalizer);

        }; // _SubfieldFactory

        const char* _SubfieldFactory::componentsVector[3] = { "_x", "_y", "_z" };
        const char* _SubfieldFactory::componentsTensor[6] = { "_xx", "_yy" "_zz", "_xy", "_yz", "_xz"};

        const char* _SubfieldFactory::genericComponent = "pylith::materials::SubfieldFactory";

    } // bc
} // pylith

// Auxiliary subfields
const std::string pylith::materials::SubfieldFactory::density = "density";
const std::string pylith::materials::SubfieldFactory::body_force = "body_force";
const std::string pylith::materials::SubfieldFactory::gravity_field = "gravitational_acceleration";
const std::string pylith::materials::SubfieldFactory::shear_modulus = "shear_modulus";
const std::string pylith::materials::SubfieldFactory::bulk_modulus = "bulk_modulus";
const std::string pylith::materials::SubfieldFactory::reference_stress = "reference_stress";
const std::string pylith::materials::SubfieldFactory::reference_strain = "reference_strain";

const std::string pylith::materials::SubfieldFactory::maxwell_time = "maxwell_time";
const std::string pylith::materials::SubfieldFactory::genmaxwell_maxwell_time = "maxwell_time";
const std::string pylith::materials::SubfieldFactory::genmaxwell_shear_modulus_ratio = "genmaxwell_shear_modulus_ratio";
const std::string pylith::materials::SubfieldFactory::powerlaw_reference_strain_rate = "powerlaw_reference_strain_rate";
const std::string pylith::materials::SubfieldFactory::powerlaw_reference_stress = "powerlaw_reference_stress";
const std::string pylith::materials::SubfieldFactory::powerlaw_exponent = "powerlaw_exponent";
const std::string pylith::materials::SubfieldFactory::total_strain = "total_strain";
const std::string pylith::materials::SubfieldFactory::deviatoric_stress = "deviatoric_sress";
const std::string pylith::materials::SubfieldFactory::viscous_strain = "viscous_strain";
const std::string pylith::materials::SubfieldFactory::genmaxwell_viscous_strain = "genmaxwell_viscous_strain";

const std::string pylith::materials::SubfieldFactory::isotropic_permeability = "isotropic_permeability";
const std::string pylith::materials::SubfieldFactory::tensor_permeability = "tensor_permeability";
const std::string pylith::materials::SubfieldFactory::drained_bulk_modulus = "drained_bulk_modulus";
const std::string pylith::materials::SubfieldFactory::biot_coefficient = "biot_coefficient";
const std::string pylith::materials::SubfieldFactory::biot_modulus = "biot_modulus";
const std::string pylith::materials::SubfieldFactory::porosity = "porosity";
const std::string pylith::materials::SubfieldFactory::solid_density = "solid_density";
const std::string pylith::materials::SubfieldFactory::fluid_density = "fluid_density";
const std::string pylith::materials::SubfieldFactory::fluid_viscosity = "fluid_viscosity";
const std::string pylith::materials::SubfieldFactory::source_density = "source_density";

// Derived subfields
const std::string pylith::materials::SubfieldFactory::cauchy_stress = "cauchy_stress";
const std::string pylith::materials::SubfieldFactory::cauchy_strain = "cauchy_strain";
const std::string pylith::materials::SubfieldFactory::bulk_density = "bulk_density";
const std::string pylith::materials::SubfieldFactory::water_content = "water_content";

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::materials::SubfieldFactory::SubfieldFactory(void) {
    GenericComponent::setName("pylith::materials::SubfieldFactory");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::materials::SubfieldFactory::~SubfieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Get subfield description.
pylith::topology::FieldBase::Description
pylith::materials::SubfieldFactory::_getDescription(const std::string& subfieldName,
                                                    const size_t spaceDim) const {
    if (subfieldName == density) {
        return _SubfieldFactory::getDensity(*_normalizer);

    } else if (subfieldName == body_force) {
        return _SubfieldFactory::getBodyForce(*_normalizer, spaceDim);

    } else if (subfieldName == gravity_field) {
        return _SubfieldFactory::getGravityField(*_normalizer, spaceDim);

    } else if (subfieldName == shear_modulus) {
        return _SubfieldFactory::getShearModulus(*_normalizer);

    } else if (subfieldName == bulk_modulus) {
        return _SubfieldFactory::getBulkModulus(*_normalizer);

    } else if (subfieldName == reference_stress) {
        return _SubfieldFactory::getReferenceStress(*_normalizer, spaceDim);

    } else if (subfieldName == reference_strain) {
        return _SubfieldFactory::getReferenceStrain(*_normalizer, spaceDim);

    } else if (subfieldName == maxwell_time) {
        return _SubfieldFactory::getMaxwellTime(*_normalizer);

    } else if (subfieldName == genmaxwell_maxwell_time) {
        return _SubfieldFactory::getGenMaxwellMaxwellTime(*_normalizer);

    } else if (subfieldName == genmaxwell_shear_modulus_ratio) {
        return _SubfieldFactory::getGenMaxwellShearModulusRatio(*_normalizer);

    } else if (subfieldName == powerlaw_reference_stress) {
        return _SubfieldFactory::getPowerlawReferenceStress(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getPowerlawExponent(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getTotalStrain(*_normalizer, spaceDim);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getDeviatoricStress(*_normalizer, spaceDim);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getViscousStrain(*_normalizer, spaceDim);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getGenMaxwellViscousStrain(*_normalizer, spaceDim);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getBiotCoefficient(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getBiotModulus(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getPorosity(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getSolidDensity(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getFluidDensity(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getFluidViscosity(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getSourceDensity(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getCauchyStress(*_normalizer, spaceDim);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getCauchyStrain(*_normalizer, spaceDim);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getPowerlawExponent(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getBulkDensity(*_normalizer);

    } else if (subfieldName == powerlaw_exponent) {
        return _SubfieldFactory::getWateContent(*_normalizer);

    } else {
        PYLITH_JOURNAL_LOGICERROR("Unrecognized subfield " << subfieldName << ".");
    } // if/else

    // Should never get here
    return pylith::topology::FieldBase::Description();
} // _getDescription


// ------------------------------------------------------------------------------------------------
// Set names of vector field components.
void
pylith::materials::_SubfieldFactory::setVectorComponentNames(pylith::topology::FieldBase::Description* description) {
    const char* componentNames[3] = {
        "_x",
        "_y",
        "_z",
    };

    const size_t numComponents = description->numComponents;
    for (int i = 0; i < numComponents; ++i) {
        description->componentNames[i] = description->label + std::string(componentNames[i]);
    } // for
} // setVectorComponentNames


// ------------------------------------------------------------------------------------------------
// Set names of vector field components.
void
pylith::materials::_SubfieldFactory::setTensorComponentNames(pylith::topology::FieldBase::Description* description) {
    const char* componentNames[6] = {
        "_xx",
        "_yy",
        "_zz",
        "_xy",
        "_yz",
        "_xz",
    };

    const size_t numComponents = description->numComponents;
    for (int i = 0; i < numComponents; ++i) {
        description->componentNames[i] = description->label + std::string(componentNames[i]);
    } // for
} // setTensorComponentNames


// ------------------------------------------------------------------------------------------------
// Get density subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getDensity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::density;
    const pylith::real densityScale = normalizer.getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getDensity


// ------------------------------------------------------------------------------------------------
// Get body force subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getBodyForce(const spatialdata::units::Nondimensional& normalizer,
                                                  const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::body_force;
    const pylith::real bodyForceScale = normalizer.getPressureScale() / normalizer.getLengthScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = spaceDim;
    description.scale = bodyForceScale;
    description.validator = nullptr;

    setVectorComponentNames(&description);

    return description;
} // getBodyForce


// ------------------------------------------------------------------------------------------------
// Get gravitational acceleration subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getGravityField(const spatialdata::units::Nondimensional& normalizer,
                                                     const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::gravity_field;

    const pylith::real lengthScale = normalizer.getLengthScale();
    const pylith::real timeScale = normalizer.getTimeScale();
    const pylith::real accelerationScale = lengthScale / (timeScale * timeScale);

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = spaceDim;
    description.scale = accelerationScale;
    description.validator = nullptr;

    setVectorComponentNames(&description);

    return description;
} // getGravityField


// ------------------------------------------------------------------------------------------------
// Get shear modulus subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getShearModulus(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::shear_modulus;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::Validators::nonnegative;

    return description;
} // getShearModulus


// ------------------------------------------------------------------------------------------------
// Get bulk modulus subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getBulkModulus(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::bulk_modulus;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getBulkModulus


// ------------------------------------------------------------------------------------------------
// Get reference stress subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getReferenceStress(const spatialdata::units::Nondimensional& normalizer,
                                                        const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::reference_stress;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = spaceDim == 2 ? 4 : 6;
    description.scale = pressureScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getReferenceStress


// ------------------------------------------------------------------------------------------------
// Get reference strain subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getReferenceStrain(const spatialdata::units::Nondimensional& normalizer,
                                                        const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::reference_strain;
    const pylith::real strainScale = 1.0;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = spaceDim == 2 ? 4 : 6;
    description.scale = strainScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getReferenceStrain


// ------------------------------------------------------------------------------------------------
// Get Maxwell time subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getMaxwellTime(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::gravity_field;
    const pylith::real timeScale = normalizer.getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    setVectorComponentNames(&description);

    return description;
} // getMaxwellTime


// ------------------------------------------------------------------------------------------------
// Get generalize Maxwell model Maxwell time subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getGenMaxwellMaxwellTime(const spatialdata::units::Nondimensional& normalizer) {
    const size_t numComponents = 3;
    const char* componentNames[3] = {
        "maxwell_time_1",
        "maxwell_time_2",
        "maxwell_time_3"
    };

    const std::string& subfieldName = pylith::materials::SubfieldFactory::genmaxwell_maxwell_time;
    const pylith::real timeScale = normalizer.getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = numComponents;
    description.componentNames.resize(numComponents);
    for (size_t i = 0; i < numComponents; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = timeScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getGenMaxwellMaxwellTime


// ------------------------------------------------------------------------------------------------
// Get generalized Maxwell model shear modulus ratio subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getGenMaxwellShearModulusRatio(const spatialdata::units::Nondimensional& normalizer) {
    const size_t numComponents = 3;
    const char* componentNames[3] = {
        "maxwell_time_1",
        "maxwell_time_2",
        "maxwell_time_3"
    };

    const std::string& subfieldName = pylith::materials::SubfieldFactory::genmaxwell_shear_modulus_ratio;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = numComponents;
    description.componentNames.resize(numComponents);
    for (size_t i = 0; i < numComponents; ++i) {
        description.componentNames[i] = componentNames[i];
    } // for
    description.scale = 1.0;
    description.validator = nullptr;

    return description;
} // getGenMaxwellShearModulusRatio


// ------------------------------------------------------------------------------------------------
// Get power law reference strain rate subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getPowerlawReferenceStrainRate(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::powerlaw_reference_strain_rate;
    const pylith::real strainRateScale = 1.0 / normalizer.getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = strainRateScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getPowerlawReferenceStrainRate


// ------------------------------------------------------------------------------------------------
// Get power law reference stress subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getPowerlawReferenceStress(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::powerlaw_reference_stress;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getPowerlawReferenceStress


// ------------------------------------------------------------------------------------------------
// Get powerlaw exponent subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getPowerlawExponent(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::powerlaw_exponent;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getPowerlawExponent


// ------------------------------------------------------------------------------------------------
// Get total strain subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getTotalStrain(const spatialdata::units::Nondimensional& normalizer,
                                                    const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::total_strain;
    const pylith::real strainScale = 1.0;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = spaceDim == 2 ? 4 : 6;
    description.scale = strainScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getTotalStrain


// ------------------------------------------------------------------------------------------------
// Get deviatoric stress subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getDeviatoricStress(const spatialdata::units::Nondimensional& normalizer,
                                                         const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::deviatoric_stress;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = spaceDim == 2 ? 4 : 6;
    description.scale = pressureScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getGravityField


// ------------------------------------------------------------------------------------------------
// Get viscous strain subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getViscousStrain(const spatialdata::units::Nondimensional& normalizer,
                                                      const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::viscous_strain;
    const pylith::real strainScale = 1.0;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = spaceDim == 2 ? 4 : 6;
    description.scale = strainScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getViscousStrain


// ------------------------------------------------------------------------------------------------
// Get generalized Maxwell model viscous strain subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getGenMaxwellViscousStrain(const spatialdata::units::Nondimensional& normalizer,
                                                                const size_t spaceDim) {
    const size_t numModels = 3;
    const char* modelNames[3] = {
        "_1",
        "_2",
        "_3",
    };
    const char* componentNames[6] = {
        "_xx",
        "_yy",
        "_zz",
        "_xy",
        "_yz",
        "_xz",
    };

    const std::string& subfieldName = pylith::materials::SubfieldFactory::genmaxwell_viscous_strain;
    const pylith::real strainScale = 1.0;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::OTHER;
    description.numComponents = 2 == spaceDim ? 4*3 : 6*3;
    description.componentNames.resize(description.numComponents);
    for (size_t iModel = 0, index = 0; iModel < numModels; ++iModel) {
        for (size_t iComponent = 0; iComponent < 6; ++iComponent) {
            description.componentNames[index++] = subfieldName + std::string(modelNames[iModel]) + std::string(componentNames[iComponent]);
        } // for
    } // for
    description.scale = strainScale;
    description.validator = nullptr;

    setVectorComponentNames(&description);

    return description;
} // getGenMaxwellViscousStrain


// ------------------------------------------------------------------------------------------------
// Get isotropic permeability subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getIsotropicPermeability(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::reference_stress;
    const pylith::real lengthScale = normalizer.getLengthScale();
    const pylith::real permeabilityScale = lengthScale*lengthScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = permeabilityScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getIsotropicPermeability


// ------------------------------------------------------------------------------------------------
// Get tensor permeability subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getTensorPermeability(const spatialdata::units::Nondimensional& normalizer,
                                                           const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::reference_stress;
    const pylith::real lengthScale = normalizer.getLengthScale();
    const pylith::real permeabilityScale = lengthScale*lengthScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = spaceDim == 2 ? 4 : 6;
    description.scale = permeabilityScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getTensorPermeability


// ------------------------------------------------------------------------------------------------
// Get drained bulk modulus subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getDrainedBulkModulus(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::reference_stress;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;
    description.validatorTolerance = 100.0;

    setTensorComponentNames(&description);

    return description;
} // getDrainedBulkModulus


// ------------------------------------------------------------------------------------------------
// Get Biot coefficient subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getBiotCoefficient(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::biot_coefficient;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = pylith::topology::FieldQuery::Validators::positive;
    description.validatorTolerance = 10.0;

    return description;
} // getBiotCoefficient


// ------------------------------------------------------------------------------------------------
// Get Biot modulus subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getBiotModulus(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::biot_coefficient;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = pressureScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;
    description.validatorTolerance = 100.0;

    return description;
} // getBiotModulus


// ------------------------------------------------------------------------------------------------
// Get gravity subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getPorosity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::porosity;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = pylith::topology::FieldQuery::Validators::nonnegative;

    return description;
} // getPorosity


// ------------------------------------------------------------------------------------------------
// Get solid density subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getSolidDensity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::solid_density;
    const pylith::real densityScale = normalizer.getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getGravityField


// ------------------------------------------------------------------------------------------------
// Get fluid density subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getFluidDensity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::fluid_density;
    const pylith::real densityScale = normalizer.getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getFluidDensity


// ------------------------------------------------------------------------------------------------
// Get fluid viscosity subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getFluidViscosity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::fluid_viscosity;
    const pylith::real pressureScale = normalizer.getPressureScale();
    const pylith::real timeScale = normalizer.getTimeScale();
    const pylith::real viscosityScale = pressureScale * timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = viscosityScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getFluidViscosity


// ------------------------------------------------------------------------------------------------
// Get source density subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getSourceDensity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::source_density;
    const pylith::real lengthScale = normalizer.getLengthScale();
    const pylith::real timeScale = normalizer.getTimeScale();
    const pylith::real sourceDensityScale = lengthScale / timeScale;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = sourceDensityScale;
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getSourceDensity


// ------------------------------------------------------------------------------------------------
// Get Cauchy stress subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getCauchyStress(const spatialdata::units::Nondimensional& normalizer,
                                                     const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::cauchy_stress;
    const pylith::real pressureScale = normalizer.getPressureScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = 2 == spaceDim ? 4 : 6;
    description.scale = pressureScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getCauchyStress


// ------------------------------------------------------------------------------------------------
// Get Cauchy strain subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getCauchyStrain(const spatialdata::units::Nondimensional& normalizer,
                                                     const size_t spaceDim) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::cauchy_strain;
    const pylith::real strainScale = 1.0;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = (3 == spaceDim) ? pylith::topology::Field::TENSOR : pylith::topology::Field::OTHER;
    description.numComponents = 2 == spaceDim ? 4 : 6;
    description.scale = strainScale;
    description.validator = nullptr;

    setTensorComponentNames(&description);

    return description;
} // getCauchyStrain


// ------------------------------------------------------------------------------------------------
// Get bulk density subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getBulkDensity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::bulk_density;
    const pylith::real densityScale = normalizer.getDensityScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = densityScale;
    description.validator = pylith::topology::FieldQuery::Validators::nonnegative;

    return description;
} // getBulkDensity


// ------------------------------------------------------------------------------------------------
// Get water content subfield description.
pylith::topology::FieldBase::Description
pylith::materials::_SubfieldFactory::getWateContent(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::materials::SubfieldFactory::water_content;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = nullptr;

    return description;
} // getWateContent


// End of file
