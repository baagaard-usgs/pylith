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
#include "TerzaghiPVE3D.hh"
#include "PressureGradientPVE3D.hh"

// Poroviscoelasticity Simple Shear Problem
TEST_CASE("SimpleShearPVE3D::TetP2P1P1::testDiscretization", "[SimpleShearPVE3D][TetP2P1P1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetP2P1P1()).testDiscretization();
}
TEST_CASE("SimpleShearPVE3D::TetP2P1P1::testResidual", "[SimpleShearPVE3D][TetP2P1P1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetP2P1P1()).testResidual();
}
TEST_CASE("SimpleShearPVE3D::TetP2P1P1::testJacobianTaylorSeries", "[SimpleShearPVE3D][TetP2P1P1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetP2P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("SimpleShearPVE3D::TetP2P1P1::testJacobianFiniteDiff", "[SimpleShearPVE3D][TetP2P1P1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::SimpleShearPVE3D::TetP2P1P1()).testJacobianFiniteDiff();
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

// Poroviscoelasticity Terzaghi problem
TEST_CASE("TerzaghiPVE3D::TetP2P1P1::testDiscretization", "[TerzaghiPVE3D][TetP2P1P1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::TetP2P1P1()).testDiscretization();
}
TEST_CASE("TerzaghiPVE3D::TetP2P1P1::testResidual", "[TerzaghiPVE3D][TetP2P1P1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::TetP2P1P1()).testResidual();
}
TEST_CASE("TerzaghiPVE3D::TetP2P1P1::testJacobianTaylorSeries", "[TerzaghiPVE3D][TetP2P1P1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::TetP2P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("TerzaghiPVE3D::TetP2P1P1::testJacobianFiniteDiff", "[TerzaghiPVE3D][TetP2P1P1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::TetP2P1P1()).testJacobianFiniteDiff();
}

TEST_CASE("TerzaghiPVE3D::HexQ2Q1Q1::testDiscretization", "[TerzaghiPVE3D][HexQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::HexQ2Q1Q1()).testDiscretization();
}
TEST_CASE("TerzaghiPVE3D::HexQ2Q1Q1::testResidual", "[TerzaghiPVE3D][HexQ2Q1Q1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::HexQ2Q1Q1()).testResidual();
}
TEST_CASE("TerzaghiPVE3D::HexQ2Q1Q1::testJacobianTaylorSeries", "[TerzaghiPVE3D][HexQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::HexQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("TerzaghiPVE3D::HexQ2Q1Q1::testJacobianFiniteDiff", "[TerzaghiPVE3D][HexQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE3D::HexQ2Q1Q1()).testJacobianFiniteDiff();
}

//Poroviscoelasticity Pressure gradient problem
TEST_CASE("PressureGradientPVE3D::TetP2P1P1::testDiscretization", "[PressureGradientPVE3D][TetP2P1P1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::TetP2P1P1()).testDiscretization();
}
TEST_CASE("PressureGradientPVE3D::TetP2P1P1::testResidual", "[PressureGradientPVE3D][TetP2P1P1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::TetP2P1P1()).testResidual();
}
TEST_CASE("PressureGradientPVE3D::TetP2P1P1::testJacobianTaylorSeries", "[PressureGradientPVE3D][TetP2P1P1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::TetP2P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientPVE3D::TetP2P1P1::testJacobianFiniteDiff", "[PressureGradientPVE3D][TetP2P1P1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::TetP2P1P1()).testJacobianFiniteDiff();
}

TEST_CASE("PressureGradientPVE3D::HexQ2Q1Q1::testDiscretization", "[PressureGradientPVE3D][HexQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::HexQ2Q1Q1()).testDiscretization();
}
TEST_CASE("PressureGradientPVE3D::HexQ2Q1Q1::testResidual", "[PressureGradientPVE3D][HexQ2Q1Q1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::HexQ2Q1Q1()).testResidual();
}
TEST_CASE("PressureGradientPVE3D::HexQ2Q1Q1::testJacobianTaylorSeries", "[PressureGradientPVE3D][HexQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::HexQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientPVE3D::HexQ2Q1Q1::testJacobianFiniteDiff", "[PressureGradientPVE3D][HexQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE3D::HexQ2Q1Q1()).testJacobianFiniteDiff();
}

// End of file
