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

#include "TestLinearPoroviscoelasticity.hh" // USES TestLinearPoroelasticity_Data

namespace pylith {
    class SimpleShearPVE3D;
}

class pylith::SimpleShearPVE3D {
public:

    // Data factory methods

    static TestLinearPoroviscoelasticity_Data* TetP2P1P1(void);

    static TestLinearPoroviscoelasticity_Data* HexQ2Q1Q1(void);

private:

    SimpleShearPVE3D(void); ///< Not implemented
}; // SimpleShear3D

// End of file