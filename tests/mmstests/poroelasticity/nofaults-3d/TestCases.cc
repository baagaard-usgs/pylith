// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** Test cases for TestLinearPoroelasticity
 */

#include "TestLinearPoroviscoelasticity.hh" // USES TestLinearPoroviscoelasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "SimpleShearPVE3D.hh"

// Poroviscoelasticity Simple Shear Problem
TEST_CASE("SimpleShearPVE3D::TetQ2Q1Q1::testDiscretization", "[SimpleShearPVE3D][TetQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetQ2Q1Q1()).testDiscretization();
}
TEST_CASE("SimpleShearPVE3D::TetQ2Q1Q1::testResidual", "[SimpleShearPVE3D][TetQ2Q1Q1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetQ2Q1Q1()).testResidual();
}
TEST_CASE("SimpleShearPVE3D::TetQ2Q1Q1::testJacobianTaylorSeries", "[SimpleShearPVE3D][TetQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("SimpleShearPVE3D::TetQ2Q1Q1::testJacobianFiniteDiff", "[SimpleShearPVE3D][TetQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetQ2Q1Q1()).testJacobianFiniteDiff();
}

TEST_CASE("SimpleShearPVE3D::HexQ2Q1Q1::testDiscretization", "[SimpleShearPVE3D][HexQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::HexQ2Q1Q1()).testDiscretization();
}
TEST_CASE("SimpleShearPVE3D::HexQ2Q1Q1::testResidual", "[SimpleShearPVE3D][HexQ2Q1Q1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::HexQ2Q1Q1()).testResidual();
}
TEST_CASE("SimpleShearPVE3D::HexQ2Q1Q1::testJacobianTaylorSeries", "[SimpleShearPVE3D][HexQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::HexQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("SimpleShearPVE3D::HexQ2Q1Q1::testJacobianFiniteDiff", "[SimpleShearPVE3D][HexQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::HexQ2Q1Q1()).testJacobianFiniteDiff();
}

// End of file
