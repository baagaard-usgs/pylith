// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
//

#include <portinfo>

#include "TestFieldMesh.hh" // Implementation of class methods

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder

#include "catch2/catch_test_macros.hpp"

namespace pylith {
    namespace topology {
        class TestFieldMesh_Cases;
    }
}

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestFieldMesh_Cases {
public:

    // Data factory methods
    static TestFieldMesh_Data* Quad(void);

};

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestFieldMesh::Quad::testConstructor", "[TestFieldMesh][Quad][testConstructor]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testConstructor();
}
TEST_CASE("TestFieldMesh::Quad::testCopyConstructor", "[TestFieldMesh][Quad][testCopyConstructor]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testCopyConstructor();
}
TEST_CASE("TestFieldMesh::Quad::testMesh", "[TestFieldMesh][Quad][testMesh]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testMesh();
}
TEST_CASE("TestFieldMesh::Quad::testGeneralAccessors", "[TestFieldMesh][Quad][testGeneralAccessors]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testGeneralAccessors();
}
TEST_CASE("TestFieldMesh::Quad::testSectionAccessors", "[TestFieldMesh][Quad][testSectionAccessors]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testSectionAccessors();
}
TEST_CASE("TestFieldMesh::Quad::testVectorAccessors", "[TestFieldMesh][Quad][testVectorAccessors]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testVectorAccessors();
}
TEST_CASE("TestFieldMesh::Quad::testSubfieldAccessors", "[TestFieldMesh][Quad][testSubfieldAccessors]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testSubfieldAccessors();
}
TEST_CASE("TestFieldMesh::Quad::testAllocate", "[TestFieldMesh][Quad][testAllocate]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testAllocate();
}
TEST_CASE("TestFieldMesh::Quad::testZeroLocal", "[TestFieldMesh][Quad][testZeroLocal]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testZeroLocal();
}
TEST_CASE("TestFieldMesh::Quad::testView", "[TestFieldMesh][Quad][testView]") {
    pylith::topology::TestFieldMesh(pylith::topology::TestFieldMesh_Cases::Quad()).testView();
}

// ------------------------------------------------------------------------------------------------
pylith::topology::TestFieldMesh_Data*
pylith::topology::TestFieldMesh_Cases::Quad(void) {
    TestFieldMesh_Data* data = new TestFieldMesh_Data();assert(data);

    const size_t numVertices = 4;
    const size_t spaceDim = 2;
    const size_t cellDim = 2;
    const size_t numCells = 1;
    const size_t numCorners = 4;
    const pylith::meshio::MeshBuilder::shape_t cellShape = pylith::meshio::MeshBuilder::QUADRILATERAL;

    static const PylithScalar vertices[numVertices*spaceDim] = {
        0.0, 0.0,
        1.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
    };
    delete data->geometry;data->geometry = new pylith::meshio::MeshBuilder::Geometry(numVertices, spaceDim, vertices);

    static const PylithInt cells[numCells*numCorners] = {
        0, 1, 3, 2,
    };
    delete data->topology;data->topology = new pylith::meshio::MeshBuilder::Topology(cellDim, numCells, numCorners, cellShape, cells);

    // Subfield A
    data->descriptionA.label = "displacement";
    data->descriptionA.vectorFieldType = FieldBase::VECTOR;
    data->descriptionA.scale = 2.0;
    data->descriptionA.numComponents = 2;
    data->descriptionA.componentNames.resize(2);
    data->descriptionA.componentNames[0] = "displacement_x";
    data->descriptionA.componentNames[1] = "displacement_y";
    data->descriptionA.validator = NULL;

    data->discretizationA.basisOrder = 1;
    data->discretizationA.quadOrder = 1;
    data->discretizationA.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;
    data->discretizationA.isBasisContinuous = true;

    static const PylithScalar _subfieldAValues[4*2] = {
        1.1, 1.2,
        2.1, 2.2,
        3.1, 3.2,
        4.1, 4.2,
    };
    data->subfieldAValues = const_cast<PylithScalar*>(_subfieldAValues);

    data->bcALabel = "boundary A";
    data->bcALabelId = 1;
    data->bcANumConstrainedDOF = 1;
    static const int _bcAConstrainedDOF[1] = { 1 };
    data->bcAConstrainedDOF = const_cast<int*>(_bcAConstrainedDOF);
    data->bcANumVertices = 2;
    static const int _bcAVertices[2] = { 1, 3, };
    data->bcAVertices = const_cast<int*>(_bcAVertices);
    data->bcAFaceValuesSize = 3;
    static const int _bcAFaceValues[3] = { 0, 1, 3 };
    data->bcAFaceValues = const_cast<int*>(_bcAFaceValues);

    // Subfield B
    data->descriptionB.label = "fluid_pressure";
    data->descriptionB.vectorFieldType = FieldBase::SCALAR;
    data->descriptionB.scale = 0.1;
    data->descriptionB.numComponents = 1;
    data->descriptionB.componentNames.resize(1);
    data->descriptionB.componentNames[0] = "fluid_pressure";
    data->descriptionB.validator = NULL;

    data->discretizationB.basisOrder = 1;
    data->discretizationB.quadOrder = 1;
    data->discretizationB.feSpace = pylith::topology::FieldBase::POLYNOMIAL_SPACE;
    data->discretizationB.isBasisContinuous = true;

    static const PylithScalar _subfieldBValues[4*1] = {
        1.3,
        2.3,
        3.3,
        4.3,
    };
    data->subfieldBValues = const_cast<PylithScalar*>(_subfieldBValues);

    data->bcBLabel = "boundary B";
    data->bcBLabelId = 1;
    data->bcBNumConstrainedDOF = 1;
    static const int _bcBConstrainedDOF[1] = { 0 };
    data->bcBConstrainedDOF = const_cast<int*>(_bcBConstrainedDOF);
    data->bcBNumVertices = 2;
    static const int _bcBVertices[2] = { 2, 3, };
    data->bcBVertices = const_cast<int*>(_bcBVertices);
    data->bcBFaceValuesSize = 3;
    static const int _bcBFaceValues[3] = { 0, 2, 3 };
    data->bcBFaceValues = const_cast<int*>(_bcBFaceValues);

    return data;
} // Quad


// End of file
