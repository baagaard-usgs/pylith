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

#include "TestSubmesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createLowerDimMesh()
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "spatialdata/geocoords/CSCart.hh" // USES CSCart

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::topology::TestSubmesh::TestSubmesh(TestSubmesh_Data* data) :
    _data(data) {
    PYLITH_METHOD_BEGIN;
    assert(data);

    _domainMesh = nullptr;
    _testMesh = nullptr;

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::TestSubmesh::~TestSubmesh(void) {
    PYLITH_METHOD_BEGIN;

    delete _domainMesh;_domainMesh = nullptr;
    delete _testMesh;_testMesh = nullptr;
    delete _data;_data = nullptr;

    PYLITH_METHOD_END;
} // destructor


// ------------------------------------------------------------------------------------------------
// Test getCoordSys(), debug(), comm().
void
pylith::topology::TestSubmesh::testAccessors(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->topology);
    assert(_data->geometry);

    _buildMesh();
    const int labelValue = 1;
    delete _testMesh;_testMesh = MeshOps::createLowerDimMesh(*_domainMesh, _data->faceGroupName, labelValue, "testAccessors");assert(_testMesh);
    assert(_testMesh->getCoordSys());

    REQUIRE(_data->topology->dimension == _testMesh->getCoordSys()->getSpaceDim());

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _testMesh->getComm(), &result);
    CHECK(int(MPI_CONGRUENT) == result);

    PYLITH_METHOD_END;
} // testAccessors


// ------------------------------------------------------------------------------------------------
// Test dimension(), numCorners(), numVertices(), numCells().
void
pylith::topology::TestSubmesh::testSizes(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->topology);
    assert(_data->geometry);

    Mesh submesh;
    CHECK(0 == submesh.getDimension());

    _buildMesh();
    const int labelValue = 1;
    delete _testMesh;_testMesh = MeshOps::createLowerDimMesh(*_domainMesh, _data->faceGroupName, labelValue, "testSizes");assert(_testMesh);
    assert(_testMesh);

    CHECK(_data->topology->dimension-1 == size_t(_testMesh->getDimension()));

    PYLITH_METHOD_END;
} // testSizes


// ------------------------------------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubmesh::testCreateLowerDimMesh(void) {
    assert(_data);
    assert(_data->topology);
    assert(_data->geometry);

    _buildMesh();
    const int labelValue = 1;
    delete _testMesh;_testMesh = MeshOps::createLowerDimMesh(*_domainMesh, _data->faceGroupName, labelValue, "testCreateLowerDimMesh");assert(_testMesh);

    REQUIRE(_data->topology->dimension-1 == size_t(_testMesh->getDimension()));

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _testMesh->getComm(), &result);
    CHECK(int(MPI_CONGRUENT) == result);

    // Check vertices
    const PetscDM dmMesh = _testMesh->getDM();assert(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const pylith::integer vStart = depthStratum.begin();
    const pylith::integer vEnd = depthStratum.end();

    const pylith::integer nvertices = _data->submeshNumVertices;
    REQUIRE(nvertices == depthStratum.size());
    for (pylith::integer v = vStart, iV = 0; v < vEnd; ++v, ++iV) {
        CHECK(_data->submeshVertices[iV] == v);
    } // for

    // Check cells
    Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
    const pylith::integer cStart = heightStratum.begin();
    const pylith::integer cEnd = heightStratum.end();

    const pylith::integer ncells = _data->submeshNumCells;
    REQUIRE(ncells == heightStratum.size());
    for (pylith::integer c = cStart, iC = 0; c < cEnd; ++c, ++iC) {
        CHECK(_data->submeshCells[iC] == c);
    } // for

    delete _testMesh;_testMesh = nullptr;
    REQUIRE_THROWS_AS(MeshOps::createLowerDimMesh(*_domainMesh, "zzyyxx", labelValue, "testFail1"), std::runtime_error);
    REQUIRE_THROWS_AS(MeshOps::createLowerDimMesh(*_domainMesh, _data->faceGroupName, labelValue+99, "testFail2"), std::runtime_error);
} // testCreateLowerDimMesh


// ------------------------------------------------------------------------------------------------
// Test constructor w/mesh.
void
pylith::topology::TestSubmesh::testCreateSubdomainMesh(void) {
    assert(_data);
    assert(_data->topology);
    assert(_data->geometry);

    _buildMesh();
    delete _testMesh;_testMesh = MeshOps::createSubdomainMesh(*_domainMesh, _data->subdomainLabel, _data->subdomainLabelValue, "Test subdomain");
    assert(_testMesh);

    REQUIRE(_data->topology->dimension == size_t(_testMesh->getDimension()));

    int result = 0;
    MPI_Comm_compare(PETSC_COMM_WORLD, _testMesh->getComm(), &result);
    CHECK(int(MPI_CONGRUENT) == result);

    // Check vertices
    const PetscDM dmMesh = _testMesh->getDM();assert(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const pylith::integer vStart = depthStratum.begin();
    const pylith::integer vEnd = depthStratum.end();

    const pylith::integer nvertices = _data->subdomainNumVertices;
    REQUIRE(nvertices == depthStratum.size());
    for (pylith::integer v = vStart, iV = 0; v < vEnd; ++v, ++iV) {
        CHECK(_data->subdomainVertices[iV] == v);
    } // for

    // Check cells
    Stratum heightStratum(dmMesh, Stratum::HEIGHT, 0);
    const pylith::integer cStart = heightStratum.begin();
    const pylith::integer cEnd = heightStratum.end();

    const pylith::integer ncells = _data->subdomainNumCells;
    REQUIRE(ncells == heightStratum.size());
    for (pylith::integer c = cStart, iC = 0; c < cEnd; ++c, ++iC) {
        CHECK(_data->subdomainCells[iC] == c);
    } // for

    delete _testMesh;_testMesh = nullptr;
    CHECK_THROWS_AS(MeshOps::createSubdomainMesh(*_domainMesh, "material-id", -9, "Test subdomain"), std::runtime_error);
    CHECK_THROWS_AS(MeshOps::createSubdomainMesh(*_domainMesh, "zzyyxx", -9, "Test subdomain"), std::runtime_error);
} // testCreateLowerDimMesh


// ------------------------------------------------------------------------------------------------
void
pylith::topology::TestSubmesh::_buildMesh(void) {
    PYLITH_METHOD_BEGIN;
    assert(_data);
    assert(_data->topology);
    assert(_data->geometry);
    assert(!_domainMesh);
    assert(!_testMesh);

    delete _domainMesh;_domainMesh = new pylith::topology::Mesh();assert(_domainMesh);
    pylith::meshio::MeshBuilder::buildMesh(_domainMesh, *_data->topology, *_data->geometry);

    spatialdata::geocoords::CSCart cs;
    cs.setSpaceDim(_data->geometry->spaceDim);
    _domainMesh->setCoordSys(&cs);

    pylith::meshio::MeshBuilder::shape_t faceShape = pylith::meshio::MeshBuilder::faceShapeFromCellShape(_data->topology->cellShape);
    pylith::integer_array faceValues(_data->faceGroup, _data->faceGroupSize);
    pylith::meshio::MeshBuilder::setFaceGroupFromCellVertices(_domainMesh, _data->faceGroupName, faceValues, faceShape);

    // Create "subdomain" by setting label of subdomain cells.
    PetscErrorCode err = PETSC_SUCCESS;
    for (size_t c = 0; c < _data->topology->numCells; ++c) {
        err = DMSetLabelValue(_domainMesh->getDM(), _data->subdomainLabel, c,
                              _data->subdomainLabelValues[c]);PYLITH_CHECK_ERROR(err);
    } // for

    PYLITH_METHOD_END;
} // _buildMesh


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::TestSubmesh_Data::TestSubmesh_Data(void) :
    topology(nullptr),
    geometry(nullptr),
    faceGroupName(nullptr),
    faceGroupSize(0),
    faceGroup(nullptr),
    submeshNumCorners(0),
    submeshNumVertices(0),
    submeshVertices(nullptr),
    submeshNumCells(0),
    submeshCells(nullptr),
    subdomainLabel(nullptr),
    subdomainLabelValues(nullptr),
    subdomainLabelValue(0),
    subdomainNumCorners(0),
    subdomainNumVertices(0),
    subdomainVertices(nullptr),
    subdomainNumCells(0),
    subdomainCells(nullptr) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::TestSubmesh_Data::~TestSubmesh_Data(void) {
    delete topology;topology = nullptr;
    delete geometry;geometry = nullptr;
}


// End of file
