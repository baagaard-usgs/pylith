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
    class TerzaghiPVE;
}

class pylith::TerzaghiPVE {
public:

    // Data factory methods

    // static TestLinearPoroviscoelasticity_Data* TriP2P1P1(void);

    // static TestLinearPoroviscoelasticity_Data* TriP3P2P2(void);

    static TestLinearPoroviscoelasticity_Data* QuadQ2Q1Q1(void);

    // static TestLinearPoroviscoelasticity_Data* QuadQ3Q2Q2(void);

private:

    TerzaghiPVE(void); ///< Not implemented
}; // PressureGradient

// End of file