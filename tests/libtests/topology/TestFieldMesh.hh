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
#include "pylith/utils/petscfwd.h" // forward declarations

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Description
#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder

namespace pylith {
    namespace topology {
        class TestFieldMesh;
        class TestFieldMesh_Data;
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestFieldMesh : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestFieldMesh(TestFieldMesh_Data* data);

    /// Destructor.
    ~TestFieldMesh(void);

    /// Test constructor.
    void testConstructor(void);

    /// Test copy constructor.
    void testCopyConstructor(void);

    /// Test mesh().
    void testMesh(void);

    /// Test getLabel(), vectorFieldType(), scale(), addDimensionOkay(), getSpaceDim().
    void testGeneralAccessors(void);

    /// Test chartSize(), getStorageSize(), selectedSection(), globalSection().
    void testSectionAccessors(void);

    /// Test localVector(), globalVector().
    void testVectorAccessors(void);

    /// Test subfieldAdd(), subfieldsSetup(), hasSubfield(), subfieldNames(), subfieldInfo().
    void testSubfieldAccessors(void);

    /// Test allocate().
    void testAllocate(void);

    /// Test zeroLocal().
    void testZeroLocal(void);

    /// Test view().
    void testView(void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    /// Initialize mesh and test field.
    void _initialize(void);

    /** Verify values in field match expected values.
     *
     * @param field Field containing values to test.
     * @param scale Scale to apply to expected values.
     */
    void _checkValues(const Field& field,
                      const PylithReal scale=1.0);

    /** Verify values in PETSc vector match expected values.
     *
     * @param vec PETSc vec containing values to test.
     * @param scale Scale to apply to expected values.
     */
    void _checkValues(const PetscVec& vec,
                      const PylithReal scale=1.0);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////
protected:

    TestFieldMesh_Data* _data; ///< Data for testing.
    Mesh* _mesh; ///< Finite-element mesh.
    Field* _field; ///< Test field associated with mesh.

}; // class TestFieldMesh

// ------------------------------------------------------------------------------------------------
class pylith::topology::TestFieldMesh_Data {
    // PUBLIC METHODS //////////////////////////////////////////////////////////
public:

    /// Constructor
    TestFieldMesh_Data(void);

    /// Destructor
    ~TestFieldMesh_Data(void);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////
public:

    pylith::meshio::MeshBuilder::Topology* topology; ///< Topology for domain mesh
    pylith::meshio::MeshBuilder::Geometry* geometry; ///< Geometry for domain mesh

    /// @defgroup Subfield A information.
    /// @{
    pylith::topology::FieldBase::Description descriptionA;
    pylith::topology::FieldBase::Discretization discretizationA;
    PylithScalar* subfieldAValues; ///< Array of values to for subfield,

    const char* bcALabel; ///< Label for boundary condition.
    int bcALabelId; ///< Label id for boundary condition.
    size_t bcANumConstrainedDOF; ///< Number of constrained DOF for boundary condition.
    int* bcAConstrainedDOF; ///< Array of constrained DOF.
    size_t bcAFaceValuesSize; ///< Number of face values (cell+vertices) associated with boundary condition.
    int* bcAFaceValues; ///< Array of face values.
    size_t bcANumVertices; ///< Number of vertices in boundary condition.
    int* bcAVertices; ///< Array of vertex indices.
    /// @}

    /// @defgroup Subfield B information.
    /// @{
    pylith::topology::FieldBase::Description descriptionB;
    pylith::topology::FieldBase::Discretization discretizationB;
    PylithScalar* subfieldBValues; ///< Array of values to for subfield,

    const char* bcBLabel; ///< Label for boundary condition.
    int bcBLabelId; ///< Label id for boundary condition.
    size_t bcBNumConstrainedDOF; ///< Number of constrained DOF for boundary condition.
    int* bcBConstrainedDOF; ///< Array of constrained DOF.
    size_t bcBFaceValuesSize; ///< Number of vertices associated with boundary condition.
    int* bcBFaceValues; ///< Array of face values.
    size_t bcBNumVertices; ///< Number of vertices in boundary condition.
    int* bcBVertices; ///< Array of vertex indices.
    /// @}

}; // TestFieldMesh_Data

// End of file
