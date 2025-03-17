// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

namespace pylith {
    namespace materials {
        // New stuff
        class Material;
        class SubfieldFactory;

        class Elasticity;
        class RheologyElasticity;
        class IsotropicLinearElasticity;
        class IsotropicLinearMaxwell;
        class IsotropicLinearGenMaxwell;
        class IsotropicPowerLaw;
        class IsotropicLinearIncompElasticity;

        class IncompressibleElasticity;
        class RheologyIncompressibleElasticity;
        class IsotropicLinearIncompElasticity;

        class Poroelasticity;
        class RheologyPoroelasticity;
        class IsotropicLinearPoroelasticity;

        class Query;

    } // materials
} // pylith

// End of file
