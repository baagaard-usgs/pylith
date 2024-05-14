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

#include "TestMeshIO.hh"

namespace pylith {
    namespace meshio {
        class TestMeshIOPetsc;
    } // meshio
} // pylith

class pylith::meshio::TestMeshIOPetsc : public TestMeshIO {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor.
    TestMeshIOPetsc(TestMeshIO_Data* data);

    /// Destructor.
    ~TestMeshIOPetsc(void);

    /// Test filename()
    void testFilename(void);

    /** Test read().
     *
     * @param[in] gmshMarkVertices If true, mark Gmsh vertices.
     */
    void testRead(const bool gmshMarkVertices=false);

    /// Test write() and read().
    void testWriteRead(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOPetsc* _io; ///< Test subject.

}; // class TestMeshIOPetsc

// End of file
