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

#include "TestFieldMesh.hh" // Implementation of class methods

#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::createDMMesh()
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor::optimizeClosure()

#include "pylith/meshio/MeshBuilder.hh" // Uses MeshBuilder
#include "pylith/utils/array.hh" // USES scalar_array

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::topology::TestFieldMesh::TestFieldMesh(TestFieldMesh_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    _mesh = nullptr;
    _field = nullptr;
    _initialize();

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::TestFieldMesh::~TestFieldMesh(void) {
    PYLITH_METHOD_BEGIN;

    delete _data;_data = nullptr;
    delete _mesh;_mesh = nullptr;
    delete _field;_field = nullptr;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldMesh::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    Mesh meshA;
    Field field(meshA);

    PYLITH_METHOD_END;
} // testConstructor


// ------------------------------------------------------------------------------------------------
// Test cloneSection().
void
pylith::topology::TestFieldMesh::testCopyConstructor(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->geometry);
    assert(_data->topology);
    assert(_mesh);
    assert(_field);

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const pylith::integer vStart = depthStratum.begin();
    const pylith::integer vEnd = depthStratum.end();

    PetscErrorCode err = PETSC_SUCCESS;

    const std::string& label = "field A";
    const std::string& fullLabel = "domain solution " + label;
    Field field(*_field);
    field.setName(label.c_str());

    const char *name = nullptr;
    err = PetscObjectGetName((PetscObject)field.getDM(), &name);assert(!err);
    CHECK(fullLabel == std::string(name));

    PetscSection section = field.getLocalSection();assert(section);
    PetscVec vec = field.getLocalVector();assert(vec);

    err = PetscObjectGetName((PetscObject) vec, &name);assert(!err);
    CHECK(fullLabel == std::string(name));

    const pylith::integer ndof = _data->descriptionA.numComponents + _data->descriptionB.numComponents;
    for (pylith::integer v = vStart, iV = 0; v < vEnd; ++v, ++iV) {
        pylith::integer dof, cdof;
        err = PetscSectionGetDof(section, v, &dof);assert(!err);
        CHECK(ndof == dof);

        err = PetscSectionGetConstraintDof(section, v, &cdof);assert(!err);

        // Count number of expected constraints on vertex.
        pylith::integer numConstraintsE = 0;
        const pylith::integer offset = _data->topology->numCells;
        for (size_t i = 0; i < _data->bcANumVertices; ++i) {
            const pylith::integer vIndex = v - offset;
            if (_data->bcAVertices[i] == vIndex) {
                numConstraintsE += _data->bcANumConstrainedDOF;
                break;
            }
        } // for
        for (size_t i = 0; i < _data->bcBNumVertices; ++i) {
            const pylith::integer vIndex = v - offset;
            if (_data->bcBVertices[i] == vIndex) {
                numConstraintsE += _data->bcBNumConstrainedDOF;
                break;
            }
        } // for
        CHECK(numConstraintsE == cdof);
    } // for

    field.deallocate();

    PYLITH_METHOD_END;
} // testCopyConstructor


// ------------------------------------------------------------------------------------------------
// Test mesh().
void
pylith::topology::TestFieldMesh::testMesh(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->topology);
    assert(_field);

    CHECK(_data->topology->dimension == size_t(_field->getMesh().getDimension()));

    PYLITH_METHOD_END;
} // testMesh


// ------------------------------------------------------------------------------------------------
// Test getLabel(), vectorFieldType(), scale(), addDimensionOkay(), getSpaceDim().
void
pylith::topology::TestFieldMesh::testGeneralAccessors(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->topology);
    assert(_field);

    // Test getLabel()
    const std::string& label = "velocity";
    const std::string& fullLabel = "domain solution " + label;
    _field->setName(label.c_str());
    CHECK(fullLabel == std::string(_field->getLabel()));
    const char* name = nullptr;
    PetscErrorCode err = PETSC_SUCCESS;
    err = PetscObjectGetName((PetscObject)_field->getDM(), &name);assert(!err);
    CHECK(fullLabel == std::string(name));

    // Test getSpaceDim()
    CHECK(size_t(_data->topology->dimension) == _field->getSpaceDim());

    PYLITH_METHOD_END;
} // testGeneralAccessors


// ------------------------------------------------------------------------------------------------
// Test chartSize(), getStorageSize(), getLocalSection(), globalSection().
void
pylith::topology::TestFieldMesh::testSectionAccessors(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->geometry);
    assert(_field);

    assert(_field->getChartSize() > 0); // vertices + edges + faces + cells
    const pylith::integer ndof = _data->descriptionA.numComponents + _data->descriptionB.numComponents;
    CHECK(_data->geometry->numVertices*ndof == size_t(_field->getStorageSize()));

    PYLITH_METHOD_END;
} // testSectionAccessors


// ------------------------------------------------------------------------------------------------
// Test getLocalVector(), getGlobalVector().
void
pylith::topology::TestFieldMesh::testVectorAccessors(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->geometry);
    assert(_data->topology);
    assert(_field);

    PetscErrorCode err = PETSC_SUCCESS;
    const char* name = nullptr;
    pylith::integer size = 0;
    const pylith::integer ndof =
        _data->geometry->numVertices * (_data->descriptionA.numComponents + _data->descriptionB.numComponents);
    const pylith::integer ndofConstrained =
        _data->bcANumConstrainedDOF*_data->bcANumVertices + _data->bcBNumConstrainedDOF*_data->bcBNumVertices;

    const PetscVec& localVec = _field->getLocalVector();assert(localVec);
    err = PetscObjectGetName((PetscObject)localVec, &name);assert(!err);
    CHECK(std::string(_field->getLabel()) == std::string(name));
    err = VecGetSize(localVec, &size);assert(!err);
    CHECK(ndof == size);

    _field->createGlobalVector();
    const PetscVec& globalVec = _field->getGlobalVector();assert(globalVec);
    _field->scatterLocalToVector(globalVec);
    err = PetscObjectGetName((PetscObject)globalVec, &name);assert(!err);
    CHECK(std::string(_field->getLabel()) == std::string(name));
    err = VecGetSize(globalVec, &size);assert(!err);
    CHECK(ndof - ndofConstrained == size);

    _field->createOutputVector();
    _field->scatterLocalToOutput();
    const PetscVec& outputVec = _field->getOutputVector();assert(outputVec);
    err = PetscObjectGetName((PetscObject)outputVec, &name);assert(!err);
    CHECK(std::string(_field->getLabel()) == std::string(name));
    err = VecGetSize(outputVec, &size);assert(!err);
    CHECK(ndof == size);
    _checkValues(outputVec);

    PYLITH_METHOD_END;
} // testVectorAccessors


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::testSubfieldAccessors(void) {
    assert(_mesh);
    assert(_field);

    // Subfields setup via subfieldAdd() and subfieldsSetup() in _initialize().

    // Test hasSubfield().
    assert(_field->hasSubfield(_data->descriptionA.label.c_str()));
    assert(_field->hasSubfield(_data->descriptionB.label.c_str()));
    assert(!_field->hasSubfield("zyxwvut987654321"));

    // Test getSubfieldNames().
    const string_vector& names = _field->getSubfieldNames();
    CHECK(size_t(2) == names.size());
    CHECK(_data->descriptionA.label == names[0]);
    CHECK(_data->descriptionB.label == names[1]);

    { // Test getSubfieldInfo() for subfield A.
        const Field::SubfieldInfo& infoA = _field->getSubfieldInfo(_data->descriptionA.label.c_str());
        CHECK(0 == infoA.index);
        CHECK(_data->descriptionA.numComponents == infoA.description.numComponents);
        CHECK(_data->descriptionA.label == infoA.description.label);
        CHECK(_data->descriptionA.vectorFieldType == infoA.description.vectorFieldType);
        CHECK(_data->descriptionA.scale == infoA.description.scale);
        const string_vector& componentNames = infoA.description.componentNames;
        REQUIRE(_data->descriptionA.numComponents == componentNames.size());
        for (size_t i = 0; i < _data->descriptionA.numComponents; ++i) {
            CHECK(_data->descriptionA.componentNames[i] == componentNames[i]);
        } // for
        CHECK(_data->discretizationA.basisOrder == infoA.fe.basisOrder);
        CHECK(_data->discretizationA.quadOrder == infoA.fe.quadOrder);
        CHECK(_data->discretizationA.feSpace == infoA.fe.feSpace);
        CHECK(_data->discretizationA.isBasisContinuous == infoA.fe.isBasisContinuous);
    } // Test getSubfieldInfo() for subfield A.

    { // Test getSubfieldInfo() for subfield B.
        const Field::SubfieldInfo& infoB = _field->getSubfieldInfo(_data->descriptionB.label.c_str());
        CHECK(1 == infoB.index);
        CHECK(_data->descriptionB.numComponents == infoB.description.numComponents);
        CHECK(_data->descriptionB.label == infoB.description.label);
        CHECK(_data->descriptionB.vectorFieldType == infoB.description.vectorFieldType);
        CHECK(_data->descriptionB.scale == infoB.description.scale);
        const string_vector& componentNames = infoB.description.componentNames;
        REQUIRE(_data->descriptionB.numComponents == componentNames.size());
        for (size_t i = 0; i < _data->descriptionB.numComponents; ++i) {
            CHECK(_data->descriptionB.componentNames[i] == componentNames[i]);
        } // for
        CHECK(_data->discretizationB.basisOrder == infoB.fe.basisOrder);
        CHECK(_data->discretizationB.quadOrder == infoB.fe.quadOrder);
        CHECK(_data->discretizationB.feSpace == infoB.fe.feSpace);
        CHECK(_data->discretizationB.isBasisContinuous == infoB.fe.isBasisContinuous);
    } // Test getSubfieldInfo() for subfield B.

    CHECK_THROWS_AS(_field->getSubfieldInfo("aabbccdd"), std::runtime_error);
} /// testSubfieldAccessors


// ------------------------------------------------------------------------------------------------
// Test allocate().
void
pylith::topology::TestFieldMesh::testAllocate(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_field);

    _checkValues(*_field);

    PYLITH_METHOD_END;
} // testAllocate


// ------------------------------------------------------------------------------------------------
// Test zeroLocal().
void
pylith::topology::TestFieldMesh::testZeroLocal(void) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);
    assert(_field);
    _field->zeroLocal();

    _checkValues(*_field, 0.0);

    PYLITH_METHOD_END;
} // testZeroLocal


// ------------------------------------------------------------------------------------------------
// Test view().
void
pylith::topology::TestFieldMesh::testView(void) {
    PYLITH_METHOD_BEGIN;
    assert(_field);
    _field->view("Testing view");

    PYLITH_METHOD_END;
} // testView


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_initialize(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->geometry);
    assert(_data->topology);

    delete _mesh;_mesh = new Mesh;assert(_mesh);
    pylith::meshio::MeshBuilder::buildMesh(_mesh, *_data->topology, *_data->geometry);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_data->geometry->spaceDim);
    _mesh->setCoordSys(&cs);

    // Setup labels
    pylith::meshio::MeshBuilder::shape_t faceShape = pylith::meshio::MeshBuilder::faceShapeFromCellShape(_data->topology->cellShape);
    pylith::integer_array faceValuesA(_data->bcAFaceValues, _data->bcAFaceValuesSize);
    pylith::meshio::MeshBuilder::setFaceGroupFromCellVertices(_mesh, _data->bcALabel, faceValuesA, faceShape);

    pylith::integer_array faceValuesB(_data->bcBFaceValues, _data->bcBFaceValuesSize);
    pylith::meshio::MeshBuilder::setFaceGroupFromCellVertices(_mesh, _data->bcBLabel, faceValuesB, faceShape);

    // Setup field
    delete _field;_field = new Field(*_mesh);
    _field->setName("solution");
    _field->subfieldAdd(_data->descriptionA, _data->discretizationA);
    _field->subfieldAdd(_data->descriptionB, _data->discretizationB);
    _field->subfieldsSetup();
    _field->createDiscretization();
    pylith::topology::CoordsVisitor::optimizeClosure(_field->getDM());

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDMLabel labelA = nullptr, labelB = nullptr;
    err = DMGetLabel(_field->getDM(), _data->bcALabel, &labelA);assert(!err);
    err = DMGetLabel(_field->getDM(), _data->bcBLabel, &labelB);assert(!err);
    const pylith::integer numLabelValues = 1;
    pylith::integer i_field = 0;
    err = DMAddBoundary(_field->getDM(), DM_BC_ESSENTIAL, "bcA", labelA, numLabelValues, &_data->bcALabelId, i_field,
                        _data->bcANumConstrainedDOF, _data->bcAConstrainedDOF, nullptr, nullptr, nullptr, nullptr);assert(!err);
    i_field = 1;
    err = DMAddBoundary(_field->getDM(), DM_BC_ESSENTIAL, "bcB", labelB, numLabelValues, &_data->bcBLabelId, i_field,
                        _data->bcBNumConstrainedDOF, _data->bcBConstrainedDOF, nullptr, nullptr, nullptr, nullptr);assert(!err);
    // Allocate field.
    _field->allocate();

    // Populate with values.
    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const pylith::integer vStart = depthStratum.begin();
    const pylith::integer vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(*_field);
    const pylith::integer numComponents = _data->descriptionA.numComponents + _data->descriptionB.numComponents;
    pylith::scalar* fieldArray = fieldVisitor.localArray();
    for (pylith::integer v = vStart, indexA = 0, indexB = 0; v < vEnd; ++v) {
        // Set values for field A
        const pylith::integer offA = fieldVisitor.sectionOffset(v);
        CHECK(numComponents == fieldVisitor.sectionDof(v));
        for (size_t d = 0; d < _data->descriptionA.numComponents; ++d) {
            fieldArray[offA+d] = _data->subfieldAValues[indexA++];
        } // for
          // Set values for field B
        const pylith::integer offB = offA + _data->descriptionA.numComponents;
        for (size_t d = 0; d < _data->descriptionB.numComponents; ++d) {
            fieldArray[offB+d] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PYLITH_METHOD_END;
} // _initialize


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_checkValues(const Field& field,
                                              const pylith::real scale) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    // Create array of values in field from subfields.
    const int numVertices = _data->geometry->numVertices;
    const int numComponentsA = _data->descriptionA.numComponents;
    const int numComponentsB = _data->descriptionB.numComponents;
    scalar_array valuesE(numVertices * (numComponentsA + numComponentsB));
    for (int iVertex = 0, index = 0, indexA = 0, indexB = 0; iVertex < numVertices; ++iVertex) {
        for (int d = 0; d < numComponentsA; ++d) {
            valuesE[index++] = _data->subfieldAValues[indexA++];
        } // for
        for (int d = 0; d < numComponentsB; ++d) {
            valuesE[index++] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const pylith::integer vStart = depthStratum.begin();
    const pylith::integer vEnd = depthStratum.end();

    VecVisitorMesh fieldVisitor(field);
    pylith::scalar* fieldArray = fieldVisitor.localArray();
    const pylith::integer numComponents = numComponentsA + numComponentsB;
    const pylith::real tolerance = 1.0e-6;
    for (pylith::integer v = vStart, index = 0; v < vEnd; ++v) {
        REQUIRE(numComponents == fieldVisitor.sectionDof(v));
        const pylith::integer off = fieldVisitor.sectionOffset(v);

        for (pylith::integer d = 0; d < numComponents; ++d) {
            CHECK_THAT(fieldArray[off+d], Catch::Matchers::WithinAbs(valuesE[index++]*scale, tolerance));
        } // for
    } // for

    PYLITH_METHOD_END;
} // _checkValues


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestFieldMesh::_checkValues(const PetscVec& vec,
                                              const pylith::real scale) {
    PYLITH_METHOD_BEGIN;
    assert(_data);

    // Create array of values in field from subfields.
    const int numVertices = _data->geometry->numVertices;
    const int numComponentsA = _data->descriptionA.numComponents;
    const int numComponentsB = _data->descriptionB.numComponents;
    scalar_array valuesE(numVertices * (numComponentsA + numComponentsB));
    for (int iVertex = 0, index = 0, indexA = 0, indexB = 0; iVertex < numVertices; ++iVertex) {
        for (int d = 0; d < numComponentsA; ++d) {
            valuesE[index++] = _data->subfieldAValues[indexA++];
        } // for
        for (int d = 0; d < numComponentsB; ++d) {
            valuesE[index++] = _data->subfieldBValues[indexB++];
        } // for
    } // for

    PetscErrorCode err = PETSC_SUCCESS;
    pylith::integer size = 0;
    PylithScalar* vecArray = nullptr;
    err = VecGetSize(vec, &size);assert(!err);
    err = VecGetArray(vec, &vecArray);assert(!err);

    const pylith::integer sizeE = numVertices * (numComponentsA + numComponentsB);
    const pylith::real tolerance = 1.0e-6;
    REQUIRE(sizeE == size);
    for (pylith::integer i = 0; i < sizeE; ++i) {
        CHECK_THAT(vecArray[i], Catch::Matchers::WithinAbs(valuesE[i]*scale, tolerance));
    } // for
    err = VecRestoreArray(vec, &vecArray);assert(!err);

    PYLITH_METHOD_END;
} // _checkValues


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::TestFieldMesh_Data::TestFieldMesh_Data(void) :
    topology(nullptr),
    geometry(nullptr),

    subfieldAValues(nullptr),
    bcALabel(nullptr),
    bcALabelId(0),
    bcANumConstrainedDOF(0),
    bcAConstrainedDOF(nullptr),
    bcAFaceValuesSize(0),
    bcAFaceValues(nullptr),
    bcANumVertices(0),
    bcAVertices(nullptr),

    subfieldBValues(nullptr),
    bcBLabel(nullptr),
    bcBLabelId(0),
    bcBNumConstrainedDOF(0),
    bcBConstrainedDOF(nullptr),
    bcBFaceValuesSize(0),
    bcBFaceValues(nullptr),
    bcBNumVertices(0),
    bcBVertices(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::TestFieldMesh_Data::~TestFieldMesh_Data(void) {
    delete topology;topology = nullptr;
    delete geometry;geometry = nullptr;
}


// End of file
