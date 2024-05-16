// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder
#include "pylith/utils/types.hh" // USES PylithScalar

namespace pylith {
    namespace topology {
        class TestSubmesh;
        class TestSubmesh_Data;
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestSubmesh : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestSubmesh(TestSubmesh_Data* data);

    /// Destructor.
    ~TestSubmesh(void);

    /// Test MeshOps::createLowerDimMesh().
    void testCreateLowerDimMesh(void);

    /// Test MeshOps::testCreateSubdomainMesh().
    void testCreateSubdomainMesh(void);

    /// Test coordsys(), debug(), comm().
    void testAccessors(void);

    /// Test dimension(), numCorners(), numVertices(), numCells(),
    void testSizes(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    // Build lower dimension mesh.
    void _buildMesh(void);

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    TestSubmesh_Data* _data; ///< Data for testing.
    Mesh* _domainMesh; ///< Mesh holding domain mesh.
    Mesh* _testMesh; ///< Test subject.

}; // class TestSubmesh

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestSubmesh_Data {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    TestSubmesh_Data(void);

    /// Destructor
    ~TestSubmesh_Data(void);

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    // GENERAL, VALUES DEPEND ON TEST CASE

    pylith::meshio::MeshBuilder::Topology* topology; ///< Topology for domain mesh
    pylith::meshio::MeshBuilder::Geometry* geometry; ///< Geometry for domain mesh

    /// @defgroup Submesh information.
    /// @{
    const char* faceGroupName; ///< Label of face group associated with submesh.
    size_t faceGroupSize; ///< Number of vertices in submesh group.
    int* faceGroup; ///< Array of vertices in submesh group.
    int submeshNumCorners; ///< Number of vertices per cell.
    int submeshNumVertices; ///< Number of vertices in submesh.
    int* submeshVertices; ///< Vertices in submesh.
    int submeshNumCells; ///< Number of cells in submesh.
    int* submeshCells; ///< Array of vertices in cells [submeshNumCells*submeshNumCorners].
    /// @}

    /// @defgroup Subdomain information.
    /// @{
    const char* subdomainLabel; ///< Label defining subdomains.
    int* subdomainLabelValues; ///< Label values defining subdomains.
    int subdomainLabelValue; ///< Label value for target subdomain.
    int subdomainNumCorners; ///< Number of vertices per cell.
    int subdomainNumVertices; ///< Number of vertices in subdomain.
    int* subdomainVertices; ///< Vertices in subdomain.
    int subdomainNumCells; ///< Number of cells in subdomain.
    int* subdomainCells; ///< Array of vertices in cells [subdomainNumCells*subdomainNumCorners].
    /// @}

}; // TestSubmesh_Data

// End of file
