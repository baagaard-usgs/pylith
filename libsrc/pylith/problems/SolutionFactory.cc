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

#include "pylith/problems/SolutionFactory.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB

#include <typeinfo> // USES "<<typeid()
#include <cassert>

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace problems {
        class _SolutionFactory {
public:

            static const char* genericComponent; ///< Name of generic PyLith component

            /** Set names of vector field components.
             *
             * @param[inout] description Description in which to set names.
             */
            static
            void setVectorComponentNames(pylith::topology::FieldBase::Description* description);

            static
            pylith::topology::FieldBase::Description getDisplacement(const spatialdata::units::Nondimensional& normalizer,
                                                                     const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getVelocity(const spatialdata::units::Nondimensional& normalizer,
                                                                 const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getPressure(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getFluidPressure(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getFluidPressureDot(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTraceStrain(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTraceStrainDot(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getLagrangeMultiplierFault(const spatialdata::units::Nondimensional& normalizer,
                                                                                const size_t spaceDim);

            static
            pylith::topology::FieldBase::Description getTemperature(const spatialdata::units::Nondimensional& normalizer);

        }; // _SubfieldFactory

        const char* _SolutionFactory::genericComponent = "pylith::problems::SolutionFactory";

    } // bc
} // pylith

const std::string pylith::problems::SolutionFactory::displacement = "displacement";
const std::string pylith::problems::SolutionFactory::velocity = "velocity";
const std::string pylith::problems::SolutionFactory::pressure = "pressure";
const std::string pylith::problems::SolutionFactory::fluid_pressure = "fluid_pressure";
const std::string pylith::problems::SolutionFactory::fluid_pressure_dot = "fluid_pressure_dot";
const std::string pylith::problems::SolutionFactory::trace_strain = "trace_strain";
const std::string pylith::problems::SolutionFactory::trace_strain_dot = "trace_strain_dot";
const std::string pylith::problems::SolutionFactory::lagrange_multiplier_fault = "lagrange_multiplier_fault";
const std::string pylith::problems::SolutionFactory::temperature = "temperature";

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::problems::SolutionFactory::SolutionFactory(void) {
    GenericComponent::setName(_SolutionFactory::genericComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::problems::SolutionFactory::~SolutionFactory(void) {}


// ------------------------------------------------------------------------------------------------
void
pylith::problems::SolutionFactory::setValues(const std::shared_ptr<spatialdata::spatialdb::SpatialDB>& db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setValues(db="<<typeid(db).name()<<")");

    _field->zeroLocal();

    pylith::topology::FieldQuery query(_field);
    query.initializeWithDefaultQueries();
    query.open(db, _normalizer->getLengthScale());
    query.query();
    query.close();

    PYLITH_METHOD_END;
} // setValues


// ------------------------------------------------------------------------------------------------
// Get subfield description.
pylith::topology::FieldBase::Description
pylith::problems::SolutionFactory::_getDescription(const std::string& subfieldName) const {
    assert(_field);
    const size_t spaceDim = _field->getSpaceDim();

    if (subfieldName == displacement) {
        return _SolutionFactory::getDisplacement(*_normalizer, spaceDim);
    } else if (subfieldName == velocity) {
        return _SolutionFactory::getVelocity(*_normalizer, spaceDim);
    } else {
        PYLITH_JOURNAL_LOGICERROR("Unrecognized subfield " << subfieldName << ".");
    } // if/else

    // Should never get here
    return pylith::topology::FieldBase::Description();
} // _getDescription


// ------------------------------------------------------------------------------------------------
// Set names of vector components in auxiliary subfield.
void
pylith::problems::_SolutionFactory::setVectorComponentNames(pylith::topology::FieldBase::Description* description) {
    assert(description);
    static const char* componentsXYZ[3] = { "_x", "_y", "_z"};

    const size_t numComponents = description->numComponents;
    description->componentNames.resize(numComponents);
    for (int i = 0; i < numComponents; ++i) {
        description->componentNames[i] = description->label + std::string(componentsXYZ[i]);
    } // for

    PYLITH_METHOD_END;
} // setVectorComponentNames


// ------------------------------------------------------------------------------------------------
// Get description of displacement.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getDisplacement(const spatialdata::units::Nondimensional& normalizer,
                                                    const size_t spaceDim) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::displacement;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = spaceDim;
    description.scale = normalizer.getLengthScale();
    description.validator = nullptr;

    setVectorComponentNames(&description);

    return description;
} // getDisplacement


// ------------------------------------------------------------------------------------------------
// Get description of velocity.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getVelocity(const spatialdata::units::Nondimensional& normalizer,
                                                const size_t spaceDim) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::velocity;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = spaceDim;
    description.scale = normalizer.getLengthScale() / normalizer.getTimeScale();
    description.validator = nullptr;

    setVectorComponentNames(&description);

    return description;
} // getVelocity


// ------------------------------------------------------------------------------------------------
// Get description of pressure.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getPressure(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::pressure;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = normalizer.getPressureScale();
    description.validator = nullptr;

    return description;
} // getPressure


// ------------------------------------------------------------------------------------------------
// Get description of fluid pressure.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getFluidPressure(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::fluid_pressure;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = normalizer.getPressureScale();
    description.validator = nullptr;

    return description;
} // getFluidPressure


// ------------------------------------------------------------------------------------------------
// Get description of time derivative of fluid pressure.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getFluidPressureDot(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::fluid_pressure_dot;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = normalizer.getPressureScale() / normalizer.getTimeScale();
    description.validator = nullptr;

    return description;
} // getFluidPressureDot


// ------------------------------------------------------------------------------------------------
// Get description of trace strain.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getTraceStrain(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::trace_strain;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0;
    description.validator = nullptr;

    return description;
} // getTraceStrain


// ------------------------------------------------------------------------------------------------
// Get description of time derivative of trace strain.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getTraceStrainDot(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::trace_strain_dot;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = 1.0 / normalizer.getTimeScale();
    description.validator = nullptr;

    return description;
} // getTraceStrainDot


// ------------------------------------------------------------------------------------------------
// Get description of Lagrange multiplier fault.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getLagrangeMultiplierFault(const spatialdata::units::Nondimensional& normalizer,
                                                               const size_t spaceDim) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::lagrange_multiplier_fault;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::VECTOR;
    description.numComponents = spaceDim;
    description.scale = normalizer.getPressureScale();
    description.validator = nullptr;

    setVectorComponentNames(&description);

    return description;
} // getLagrangeMultiplierFault


// ------------------------------------------------------------------------------------------------
// Get description of temperature.
pylith::topology::FieldBase::Description
pylith::problems::_SolutionFactory::getTemperature(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::problems::SolutionFactory::temperature;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = normalizer.getTemperatureScale();
    description.validator = nullptr;

    return description;
} // getTemperature


// End of file
