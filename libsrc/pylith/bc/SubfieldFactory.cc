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

#include "pylith/bc/SubfieldFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class _SubfieldFactory {
public:

            ///< Names of field components in XYZ coordinate system.
            static const char* componentsXYZ[3];

            ///< Names of field components in 2D tangential/normal coordinate system.
            static const char* componentsTN[2];

            ///< Names of field components in 3D tangential/normal coordinate system.
            static const char* componentsTTN[3];

            static const char* genericComponent; ///< Name of generic PyLith component

            /** Set names of vector field components.
             *
             * @param[inout] description Description in which to set names.
             * @param[in] spaceDim Spatial dimension.
             */
            static
            void setVectorComponentNames(pylith::topology::FieldBase::Description* description,
                                         const pylith::bc::SubfieldFactory::ComponentsType componentsType);

            static
            pylith::topology::FieldBase::Description getDensity(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getVp(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getVs(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getInitialAmplitude(const pylith::topology::FieldBase::Description& refDescription,
                                                                         const pylith::bc::SubfieldFactory::ComponentsType componentsType,
                                                                         const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getRateAmplitude(const pylith::topology::FieldBase::Description& refDescription,
                                                                      const pylith::bc::SubfieldFactory::ComponentsType componentsType,
                                                                      const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getRateStartTime(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTimeHistoryAmplitude(const pylith::topology::FieldBase::Description& refDescription,
                                                                             const pylith::bc::SubfieldFactory::ComponentsType componentsType,
                                                                             const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTimeHistoryStartTime(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTimeHistoryValue(const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getNormalDir(const size_t spaceDim,
                                                                  const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTangentialDirHoriz(const size_t spaceDim,
                                                                           const spatialdata::units::Nondimensional& normalizer);

            static
            pylith::topology::FieldBase::Description getTangentialDirVert(const size_t spaceDim,
                                                                          const spatialdata::units::Nondimensional& normalizer);

        }; // _SubfieldFactory

        const char* _SubfieldFactory::componentsXYZ[3] = { "_x", "_y", "_z" };
        const char* _SubfieldFactory::componentsTN[2] = { "_tangential", "_normal" };
        const char* _SubfieldFactory::componentsTTN[3] = { "_tangential_1", "_tangential_2", "_normal" };

        const char* _SubfieldFactory::genericComponent = "pylith::bc::SubfieldFactory";

    } // bc
} // pylith

// Auxiliary subfields
const std::string pylith::bc::SubfieldFactory::density = "density";
const std::string pylith::bc::SubfieldFactory::vp = "vp";
const std::string pylith::bc::SubfieldFactory::vs = "vs";
const std::string pylith::bc::SubfieldFactory::initial_amplitude = "initial_amplitude";
const std::string pylith::bc::SubfieldFactory::rate_amplitude = "rate_amplitude";
const std::string pylith::bc::SubfieldFactory::rate_start_time = "rate_start_time";
const std::string pylith::bc::SubfieldFactory::time_history_amplitude = "time_history_amplitude";
const std::string pylith::bc::SubfieldFactory::time_history_start_time = "time_history_start_time";
const std::string pylith::bc::SubfieldFactory::time_history_value = "time_history_value";

// Diagnostic subfields
const std::string pylith::bc::SubfieldFactory::normal_dir = "normal_dir";
const std::string pylith::bc::SubfieldFactory::tangential_dir_horiz = "tangential_dir_horiz";
const std::string pylith::bc::SubfieldFactory::tangential_dir_vert = "tangential_dir_vert";

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::SubfieldFactory::SubfieldFactory(const ComponentsType componentsType) :
    _componentsType(componentsType) {
    GenericComponent::setName(_SubfieldFactory::genericComponent);
} // constructor


// ------------------------------------------------------------------------------------------------
// Set description of reference field.
void
pylith::bc::SubfieldFactory::setRefDescription(const std::shared_ptr<pylith::topology::FieldBase::Description>& refDescription) {
    _refDescription = refDescription;
} // setRefDescription


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::SubfieldFactory::~SubfieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Get subfield description.
pylith::topology::FieldBase::Description
pylith::bc::SubfieldFactory::_getDescription(const std::string& subfieldName,
                                             const size_t spaceDim) const {
    if (subfieldName == density) {
        return _SubfieldFactory::getDensity(*_normalizer);
    } else if (subfieldName == vp) {
        return _SubfieldFactory::getVp(*_normalizer);
    } else if (subfieldName == vs) {
        return _SubfieldFactory::getVs(*_normalizer);
    } else if (subfieldName == initial_amplitude) {
        return _SubfieldFactory::getInitialAmplitude(*_refDescription, _componentsType, *_normalizer);
    } else if (subfieldName == rate_amplitude) {
        return _SubfieldFactory::getRateAmplitude(*_refDescription, _componentsType, *_normalizer);
    } else if (subfieldName == rate_start_time) {
        return _SubfieldFactory::getRateStartTime(*_normalizer);
    } else if (subfieldName == time_history_amplitude) {
        return _SubfieldFactory::getTimeHistoryAmplitude(*_refDescription, _componentsType, *_normalizer);
    } else if (subfieldName == time_history_start_time) {
        return _SubfieldFactory::getTimeHistoryStartTime(*_normalizer);
    } else if (subfieldName == time_history_value) {
        return _SubfieldFactory::getTimeHistoryValue(*_normalizer);

    } else if (subfieldName == normal_dir) {
        return _SubfieldFactory::getNormalDir(spaceDim, *_normalizer);
    } else if (subfieldName == tangential_dir_horiz) {
        return _SubfieldFactory::getTangentialDirHoriz(spaceDim, *_normalizer);
    } else if (subfieldName == tangential_dir_vert) {
        return _SubfieldFactory::getTangentialDirVert(spaceDim, *_normalizer);
    } else {
        PYLITH_JOURNAL_LOGICERROR("Unrecognized subfield " << subfieldName << ".");
    } // if/else

    // Should never get here
    return pylith::topology::FieldBase::Description();
} // _getDescription


// ------------------------------------------------------------------------------------------------
// Set names of vector components in subfield.
void
pylith::bc::_SubfieldFactory::setVectorComponentNames(pylith::topology::FieldBase::Description* description,
                                                      const pylith::bc::SubfieldFactory::ComponentsType componentsType) {
    assert(description);

    const char** componentNames = nullptr;
    const size_t numComponents = description->numComponents;
    switch (componentsType) {
    case pylith::bc::SubfieldFactory::XYZ:
        componentNames = componentsXYZ;
        break;
    case pylith::bc::SubfieldFactory::TANGENTIAL_NORMAL:
        if (2 == numComponents) {
            componentNames = componentsTN;
        } else if (3 == numComponents) {
            componentNames = componentsTTN;
        } // if/else
        break;
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown case for auxiliary component reference ("<<componentsType<<") and " << numComponents << " components.");
    } // if

    for (int i = 0; i < numComponents; ++i) {
        description->componentNames[i] = description->label + std::string(componentNames[i]);
    } // for

    PYLITH_METHOD_END;
} // _setVectorComponentNames


// ---------------------------------------------------------------------------------------------------------------------
// Get description of density.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getDensity(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::density;

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = normalizer.getDensityScale();
    description.validator = pylith::topology::FieldQuery::Validators::positive;

    return description;
} // getDensity


// ---------------------------------------------------------------------------------------------------------------------
// Get description of shear wave speed.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getVs(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::vs;
    const pylith::real velocityScale = normalizer.getLengthScale() / normalizer.getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = velocityScale;
    description.validator = nullptr;

    return description;
} // getVs


// ---------------------------------------------------------------------------------------------------------------------
// Get description of dilatational wave speed.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getVp(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::vp;
    const pylith::real velocityScale = normalizer.getLengthScale() / normalizer.getTimeScale();

    pylith::topology::Field::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = velocityScale;
    description.validator = nullptr;

    return description;
} // getVp


// ------------------------------------------------------------------------------------------------
// Get description of initial amplitude.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getInitialAmplitude(const pylith::topology::FieldBase::Description& refDescription,
                                                  const pylith::bc::SubfieldFactory::ComponentsType componentsType,
                                                  const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::initial_amplitude;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = refDescription.vectorFieldType;
    description.numComponents = refDescription.numComponents;
    description.scale = refDescription.scale;
    description.validator = nullptr;

    switch (description.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        description.componentNames.resize(1);
        description.componentNames[0] = subfieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        setVectorComponentNames(&description, componentsType);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown vector field case.");
    } // switch

    return description;
} // getInitialAmplitude


// ------------------------------------------------------------------------------------------------
// Get description of rate amplitude.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getRateAmplitude(const pylith::topology::FieldBase::Description& refDescription,
                                               const pylith::bc::SubfieldFactory::ComponentsType componentsType,
                                               const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::rate_amplitude;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = refDescription.vectorFieldType;
    description.numComponents = refDescription.numComponents;
    description.scale = refDescription.scale / normalizer.getTimeScale();
    description.validator = nullptr;

    switch (description.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        description.componentNames.resize(1);
        description.componentNames[0] = subfieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        setVectorComponentNames(&description, componentsType);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown vector field case.");
    } // switch

    return description;
} // getRateAmplitude


// ------------------------------------------------------------------------------------------------
// Get description of rate start time.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getRateStartTime(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::rate_start_time;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = normalizer.getTimeScale();
    description.validator = nullptr;

    return description;
} // getRateStartTime


// ------------------------------------------------------------------------------------------------
// Get description of time history amplitude.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getTimeHistoryAmplitude(const pylith::topology::FieldBase::Description& refDescription,
                                                      const pylith::bc::SubfieldFactory::ComponentsType componentsType,
                                                      const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::rate_amplitude;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = refDescription.vectorFieldType;
    description.numComponents = refDescription.numComponents;
    description.scale = refDescription.scale;
    description.validator = nullptr;

    switch (description.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        description.componentNames.resize(1);
        description.componentNames[0] = subfieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        setVectorComponentNames(&description, componentsType);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_LOGICERROR("Unknown vector field case.");
    } // switch

    return description;
} // getTimeHistoryAmplitude


// ------------------------------------------------------------------------------------------------
// Get description of time history start time.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getTimeHistoryStartTime(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::time_history_start_time;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.scale = normalizer.getTimeScale();
    description.validator = nullptr;

    return description;
} // getTimeHistoryStartTime


// ------------------------------------------------------------------------------------------------
// Get description of time history value.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getTimeHistoryValue(const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::time_history_value;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = subfieldName;
    description.validator = nullptr;

    return description;
} // getTimeHistoryValue


// ------------------------------------------------------------------------------------------------
// Get description of normal direction.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getNormalDir(const size_t spaceDim,
                                           const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::normal_dir;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::VECTOR;
    description.numComponents = spaceDim;
    description.scale = 1.0;
    description.validator = nullptr;

    setVectorComponentNames(&description, pylith::bc::SubfieldFactory::XYZ);

    return description;
} // getNormalDir


// ------------------------------------------------------------------------------------------------
// Get description of horizontal tangential direction.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getTangentialDirHoriz(const size_t spaceDim,
                                                    const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::tangential_dir_horiz;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::VECTOR;
    description.numComponents = spaceDim;
    description.scale = 1.0;
    description.validator = nullptr;

    setVectorComponentNames(&description, pylith::bc::SubfieldFactory::XYZ);

    return description;
} // getTangentialDirHoriz


// ------------------------------------------------------------------------------------------------
// Get description of vertical tangential direction.
pylith::topology::FieldBase::Description
pylith::bc::_SubfieldFactory::getTangentialDirVert(const size_t spaceDim,
                                                   const spatialdata::units::Nondimensional& normalizer) {
    const std::string& subfieldName = pylith::bc::SubfieldFactory::tangential_dir_vert;

    pylith::topology::FieldBase::Description description;
    description.label = subfieldName;
    description.alias = subfieldName;
    description.vectorFieldType = pylith::topology::FieldBase::VECTOR;
    description.numComponents = spaceDim;
    description.scale = 1.0;
    description.validator = nullptr;

    setVectorComponentNames(&description, pylith::bc::SubfieldFactory::XYZ);

    return description;
} // getTangentialDirVert


// End of file
