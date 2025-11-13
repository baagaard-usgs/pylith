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

#include "pylith/meshio/OutputSubfield.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/FieldOps.hh" // USES FieldOps
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/RefineInterpolator.hh" // USES RefineInterpolator
#include "pylith/fekernels/Solution.hh" // USES Solution::passThruSubfield

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "petscdm.h" // USES DMReorderSectionSetDefault()

#include <typeinfo> // USES typeid()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        class _OutputSubfield {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static pylith::integer create;
                static pylith::integer createBasisOrder;
                static pylith::integer setLabel;
                static pylith::integer project;
                static pylith::integer projectWithLabel;
                static pylith::integer extractSubfield;
            };

        }; // _OutputSubfield
    } // meshio
} // pylith

pylith::utils::EventLogger pylith::meshio::_OutputSubfield::Events::logger;
pylith::integer pylith::meshio::_OutputSubfield::Events::create;
pylith::integer pylith::meshio::_OutputSubfield::Events::createBasisOrder;
pylith::integer pylith::meshio::_OutputSubfield::Events::setLabel;
pylith::integer pylith::meshio::_OutputSubfield::Events::project;
pylith::integer pylith::meshio::_OutputSubfield::Events::projectWithLabel;
pylith::integer pylith::meshio::_OutputSubfield::Events::extractSubfield;

// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_OutputSubfield::Events::init(void) {
    logger.setClassName("OutputSubfield");
    logger.initialize();
    create = logger.registerEvent("PL:OutputSubfield:create");
    createBasisOrder = logger.registerEvent("PL:OutputSubfield:createBasisOrder");
    setLabel = logger.registerEvent("PL:OutputSubfield:setLabel");
    project = logger.registerEvent("PL:OutputSubfield:project");
    projectWithLabel = logger.registerEvent("PL:OutputSubfield:projectWithLabel");
    extractSubfield = logger.registerEvent("PL:OutputSubfield:extractSubfield");
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSubfield::OutputSubfield(void) :
    _subfieldIndex(-1),
    _labelValue(0),
    _projectDM(nullptr),
    _projectVector(nullptr),
    _projectVectorInterp(nullptr),
    _fn(pylith::fekernels::Solution::passThruSubfield),
    _outputDM(nullptr),
    _outputVector(nullptr) {
    _OutputSubfield::Events::init();
}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSubfield::~OutputSubfield(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSubfield::deallocate(void) {
    _description.reset();
    _discretization.reset();
    _interpolator.reset();

    PetscErrorCode err = PETSC_SUCCESS;
    err = VecDestroy(&_projectVector);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_projectVectorInterp);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&_projectDM);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_outputVector);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&_outputDM);PYLITH_CHECK_ERROR(err);
} // deallocate


// ------------------------------------------------------------------------------------------------
// Create OutputSubfield from Field.
std::unique_ptr<pylith::meshio::OutputSubfield>
pylith::meshio::OutputSubfield::create(const pylith::topology::Field& field,
                                       const pylith::topology::Mesh& mesh,
                                       const char* name,
                                       const int basisOrder,
                                       const int refineLevels) {
    PYLITH_METHOD_BEGIN;
    // _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::create);

    std::unique_ptr<OutputSubfield> subfield(new OutputSubfield);assert(subfield);

    const pylith::topology::Field::SubfieldInfo& info = field.getSubfieldInfo(name);
    subfield->_subfieldIndex = info.index;
    subfield->_description = std::make_unique<pylith::topology::Field::Description>(info.description);
    const int outputBasisOrder = std::min(basisOrder, info.fe.basisOrder);

    // Discretization for projection
    pylith::topology::FieldBase::Discretization projectDiscretization = info.fe;
    projectDiscretization.dimension = mesh.getDimension();
    projectDiscretization.basisOrder = refineLevels ? info.fe.basisOrder : outputBasisOrder;

    // Discretization for output
    subfield->_discretization = std::make_unique<pylith::topology::FieldBase::Discretization>(projectDiscretization);
    subfield->_discretization->basisOrder = outputBasisOrder;

    PetscErrorCode err = PETSC_SUCCESS;

    // Setup PETSc DM for projection
    const char* meshName = nullptr;
    err = PetscObjectGetName((PetscObject) mesh.getDM(), &meshName);PYLITH_CHECK_ERROR(err);
    const std::string& projectName = meshName + std::string(" ") + std::string(name);
    err = DMClone(mesh.getDM(), &subfield->_projectDM);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectDM, projectName.c_str());PYLITH_CHECK_ERROR(err);
    err = DMReorderSectionSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE);PYLITH_CHECK_ERROR(err);
    err = DMReorderSectionSetType(subfield->_projectDM, nullptr);PYLITH_CHECK_ERROR(err);
    err = DMPlexReorderSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE);

    // Setup PETSc FE (discretization) for projection
    PetscFE projectFE = pylith::topology::FieldOps::createFE(projectDiscretization, subfield->_projectDM,
                                                             info.description.numComponents);assert(projectFE);
    err = PetscFESetName(projectFE, info.description.label.c_str());PYLITH_CHECK_ERROR(err);
    err = DMSetField(subfield->_projectDM, 0, nullptr, (PetscObject)projectFE);PYLITH_CHECK_ERROR(err);
    err = DMSetFieldAvoidTensor(subfield->_projectDM, 0, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
    err = PetscFEDestroy(&projectFE);PYLITH_CHECK_ERROR(err);
    err = DMCreateDS(subfield->_projectDM);PYLITH_CHECK_ERROR(err);

    if (!refineLevels) {
        subfield->_outputDM = subfield->_projectDM;
        err = PetscObjectReference((PetscObject)subfield->_outputDM);PYLITH_CHECK_ERROR(err);
    } else {
        subfield->_interpolator = std::make_unique<pylith::topology::RefineInterpolator>();
        assert(subfield->_interpolator);
        subfield->_interpolator->initialize(subfield->_projectDM, refineLevels, outputBasisOrder, info.description, *subfield->_discretization);
        subfield->_outputDM = subfield->_interpolator->getOutputDM();
        err = PetscObjectReference((PetscObject)subfield->_outputDM);PYLITH_CHECK_ERROR(err);
    } // if/else

    err = DMCreateGlobalVector(subfield->_projectDM, &subfield->_projectVector);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectVector, name);PYLITH_CHECK_ERROR(err);
    if (refineLevels) {
        err = DMCreateGlobalVector(subfield->_outputDM, &subfield->_outputVector);PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject)subfield->_outputVector, name);PYLITH_CHECK_ERROR(err);
        err = VecDuplicate(subfield->_projectVector, &subfield->_projectVectorInterp);PYLITH_CHECK_ERROR(err);
    } else {
        subfield->_outputVector = subfield->_projectVector;
        err = PetscObjectReference((PetscObject)subfield->_outputVector);PYLITH_CHECK_ERROR(err);
    }

    // _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::create);
    PYLITH_METHOD_RETURN(subfield);
}


// ------------------------------------------------------------------------------------------------
// Create OutputSubfield from Field.
std::unique_ptr<pylith::meshio::OutputSubfield>
pylith::meshio::OutputSubfield::create(const pylith::topology::Field& field,
                                       const pylith::topology::Mesh& mesh,
                                       const char* name) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::createBasisOrder);

    std::unique_ptr<OutputSubfield> subfield(new OutputSubfield());assert(subfield);

    const pylith::topology::Field::SubfieldInfo& info = field.getSubfieldInfo(name);
    subfield->_subfieldIndex = info.index;
    subfield->_description = std::make_unique<pylith::topology::Field::Description>(info.description);

    PetscErrorCode err = PETSC_SUCCESS;
    err = DMClone(mesh.getDM(), &subfield->_projectDM);PYLITH_CHECK_ERROR(err);assert(subfield->_projectDM);
    err = DMReorderSectionSetDefault(subfield->_projectDM, DM_REORDER_DEFAULT_FALSE);PYLITH_CHECK_ERROR(err);
    err = DMReorderSectionSetType(subfield->_projectDM, nullptr);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectDM, name);PYLITH_CHECK_ERROR(err);

    pylith::topology::VecVisitorMesh fieldVisitor(field, name);

    PetscSection subfieldSection = nullptr;
    pylith::integer pStart = 0, pEnd = 0;
    err = PetscSectionClone(fieldVisitor.selectedSection(), &subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetChart(fieldVisitor.selectedSection(), &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    for (pylith::integer point = pStart, offset = 0; point < pEnd; ++point) {
        const pylith::integer numDof = fieldVisitor.sectionDof(point);
        err = PetscSectionSetOffset(subfieldSection, point, offset);PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetDof(subfieldSection, point, numDof);PYLITH_CHECK_ERROR(err);
        offset += numDof;
    } // for
    err = DMSetLocalSection(subfield->_projectDM, subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&subfieldSection);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(subfield->_projectDM, &subfield->_projectVector);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject)subfield->_projectVector, name);PYLITH_CHECK_ERROR(err);

    subfield->_outputDM = subfield->_projectDM;
    err = PetscObjectReference((PetscObject) subfield->_outputDM);PYLITH_CHECK_ERROR(err);
    subfield->_outputVector = subfield->_projectVector;
    err = PetscObjectReference((PetscObject) subfield->_outputVector);PYLITH_CHECK_ERROR(err);

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::createBasisOrder);
    PYLITH_METHOD_RETURN(subfield);
}


// ------------------------------------------------------------------------------------------------
// Set label name and value.
void
pylith::meshio::OutputSubfield::setName(const char* name,
                                        const int value) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::setLabel);

    _labelName = name;
    _labelValue = value;

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::setLabel);
    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Get description of subfield.
const pylith::topology::FieldBase::Description&
pylith::meshio::OutputSubfield::getDescription(void) const {
    return *_description;
}


// ------------------------------------------------------------------------------------------------
// Get basis order of subfield.
int
pylith::meshio::OutputSubfield::getBasisOrder(void) const {
    return _discretization->basisOrder;
}


// ------------------------------------------------------------------------------------------------
// Get filtered PETSc global vector.
PetscVec
pylith::meshio::OutputSubfield::getOutputVector(void) const {
    return _outputVector;
}


// ------------------------------------------------------------------------------------------------
// Get PETSc DM for filtered vector.
PetscDM
pylith::meshio::OutputSubfield::getOutputDM(void) const {
    return _outputDM;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield data from global PETSc vector with subfields.
void
pylith::meshio::OutputSubfield::project(const PetscVec& fieldVector) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::project);
    assert(fieldVector);
    assert(_projectVector);
    assert(_outputVector);

    PetscErrorCode err = PETSC_SUCCESS;
    const pylith::real t = pylith::real(_subfieldIndex) + 0.01; // :KLUDGE: Easiest way to get subfield to extract into
                                                                // fn

    err = DMProjectField(_projectDM, t, fieldVector, &_fn, INSERT_VALUES, _projectVector);PYLITH_CHECK_ERROR(err);
    if (_interpolator) {
        pylith::topology::FieldOps::transformVector(&_projectVectorInterp, _interpolator->getInputDM(), _projectVector, _projectDM);
        _interpolator->interpolate(&_outputVector, _projectVectorInterp);
    } // if
    err = VecScale(_outputVector, _description->scale);PYLITH_CHECK_ERROR(err);

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::project);
    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield data from global PETSc vector with subfields.
void
pylith::meshio::OutputSubfield::projectWithLabel(const PetscVec& fieldVector) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::projectWithLabel);
    assert(fieldVector);
    assert(_projectVector);
    assert(_outputVector);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDMLabel dmLabel = nullptr;
    err = DMGetLabel(_projectDM, _labelName.c_str(), &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelComplete(_projectDM, dmLabel);

    const pylith::real t = pylith::real(_subfieldIndex) + 0.01; // :KLUDGE: Easiest way to get subfield to extract into
                                                                // fn
    err = DMProjectFieldLabel(_projectDM, t, dmLabel, 1, &_labelValue, PETSC_DETERMINE, nullptr, fieldVector, &_fn, INSERT_VALUES, _projectVector);PYLITH_CHECK_ERROR(err);
    if (_interpolator) {
        pylith::topology::FieldOps::transformVector(&_projectVectorInterp, _interpolator->getInputDM(), _projectVector, _projectDM);
        _interpolator->interpolate(&_outputVector, _projectVectorInterp);
    } // if
    err = VecScale(_outputVector, _description->scale);PYLITH_CHECK_ERROR(err);

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::projectWithLabel);
    PYLITH_METHOD_END;
}


// ------------------------------------------------------------------------------------------------
// Extract subfield from field.
void
pylith::meshio::OutputSubfield::extractSubfield(const pylith::topology::Field& field,
                                                const pylith::integer subfieldIndex) {
    PYLITH_METHOD_BEGIN;
    _OutputSubfield::Events::logger.eventBegin(_OutputSubfield::Events::extractSubfield);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscSection subfieldSection = nullptr;
    pylith::integer storageSize = 0;
    err = PetscSectionGetField(field.getLocalSection(), subfieldIndex, &subfieldSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionGetStorageSize(subfieldSection, &storageSize);PYLITH_CHECK_ERROR(err);

    PetscVec subfieldVector = this->getOutputVector();
    pylith::integer subfieldSize = 0;
    err = VecGetLocalSize(subfieldVector, &subfieldSize);PYLITH_CHECK_ERROR(err);
    assert(subfieldSize == storageSize);

    pylith::integer pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(subfieldSection, &pStart, &pEnd);

    pylith::topology::VecVisitorMesh fieldVisitor(field);
    pylith::scalar* solnArray = fieldVisitor.localArray();
    pylith::scalar* subfieldArray = nullptr;
    err = VecGetArray(subfieldVector, &subfieldArray);PYLITH_CHECK_ERROR(err);

    for (pylith::integer point = pStart, indexVec = 0; point < pEnd; ++point) {
        const pylith::integer solnOffset = fieldVisitor.sectionSubfieldOffset(subfieldIndex, point);
        const pylith::integer solnDof = fieldVisitor.sectionSubfieldDof(subfieldIndex, point);

        for (pylith::integer iDof = 0; iDof < solnDof; ++iDof) {
            // Dimensionalize values while extracting subfield.
            subfieldArray[indexVec++] = solnArray[solnOffset+iDof] * _description->scale;
        } // for
    } // for

    err = VecRestoreArray(subfieldVector, &subfieldArray);PYLITH_CHECK_ERROR(err);

    _OutputSubfield::Events::logger.eventEnd(_OutputSubfield::Events::extractSubfield);
    PYLITH_METHOD_END;
} // extractSubfield


// End of file
