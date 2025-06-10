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

#include "TestLinearPoroelasticity.hh" // USES TestLineaerPoroelasticity
#include "TestLinearPoroviscoelasticity.hh" // USES TestLinearPoroviscoelasticity

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
#include "PressureGradient.hh"
#include "PressureGradientPVE.hh"
#include "TerzaghiPVE.hh"

/// Poroelasticity Pressure Gradient

// TriP2P1P1
TEST_CASE("PressureGradient::TriP2P1P1::testDiscretization", "[PressureGradient][TriP2P1P1][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP2P1P1::testResidual", "[PressureGradient][TriP2P1P1][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testResidual();
}
TEST_CASE("PressureGradient::TriP2P1P1::testJacobianTaylorSeries", "[PressureGradient][TriP2P1P1][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP2P1P1::testJacobianFiniteDiff", "[PressureGradient][TriP2P1P1][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1()).testJacobianFiniteDiff();
}

// TriP3P2P2
TEST_CASE("PressureGradient::TriP3P2P2::testDiscretization", "[PressureGradient][TriP3P2P2][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP3P2P2::testResidual", "[PressureGradient][TriP3P2P2][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testResidual();
}
TEST_CASE("PressureGradient::TriP3P2P2::testJacobianTaylorSeries", "[PressureGradient][TriP3P2P2][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP3P2P2::testJacobianFiniteDiff", "[PressureGradient][TriP3P2P2][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testDiscretization", "[PressureGradient][QuadQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testResidual", "[PressureGradient][QuadQ2Q1Q1][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testJacobianTaylorSeries", "[PressureGradient][QuadQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1::testJacobianFiniteDiff", "[PressureGradient][QuadQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testDiscretization", "[PressureGradient][QuadQ3Q2Q2][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testResidual", "[PressureGradient][QuadQ3Q2Q2][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testJacobianTaylorSeries", "[PressureGradient][QuadQ3Q2Q2][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2::testJacobianFiniteDiff", "[PressureGradient][QuadQ3Q2Q2][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2()).testJacobianFiniteDiff();
}

// TriP2P1P1 w/state variables
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testDiscretization", "[PressureGradient][TriP2P1P1_StateVars][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testResidual", "[PressureGradient][TriP2P1P1_StateVars][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testJacobianTaylorSeries", "[PressureGradient][TriP2P1P1_StateVars][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP2P1P1_StateVars::testJacobianFiniteDiff", "[PressureGradient][TriP2P1P1_StateVars][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP2P1P1_StateVars()).testJacobianFiniteDiff();
}

// TriP3P2P2 with state variables
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testDiscretization", "[PressureGradient][TriP3P2P2_StateVars][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testResidual", "[PressureGradient][TriP3P2P2_StateVars][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testJacobianTaylorSeries", "[PressureGradient][TriP3P2P2_StateVars][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::TriP3P2P2_StateVars::testJacobianFiniteDiff", "[PressureGradient][TriP3P2P2_StateVars][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::TriP3P2P2_StateVars()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1 with state variables
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testDiscretization", "[PressureGradient][QuadQ2Q1Q1_StateVars][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testResidual", "[PressureGradient][QuadQ2Q1Q1_StateVars][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testJacobianTaylorSeries", "[PressureGradient][QuadQ2Q1Q1_StateVars][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ2Q1Q1_StateVars::testJacobianFiniteDiff", "[PressureGradient][QuadQ2Q1Q1_StateVars][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ2Q1Q1_StateVars()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2 with state variables
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testDiscretization", "[PressureGradient][QuadQ3Q2Q2][discretization]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testDiscretization();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testResidual", "[PressureGradient][QuadQ3Q2Q2][residual]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testResidual();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testJacobianTaylorSeries", "[PressureGradient][QuadQ3Q2Q2][Jacobian Taylor series]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradient::QuadQ3Q2Q2_StateVars::testJacobianFiniteDiff", "[PressureGradient][QuadQ3Q2Q2][Jacobian finite difference]") {
    pylith::TestLinearPoroelasticity(pylith::PressureGradient::QuadQ3Q2Q2_StateVars()).testJacobianFiniteDiff();
}

/// Poroviscoelasticity Pressure Gradient 

// TriP2P1P1
TEST_CASE("PressureGradientPVE::TriP2P1P1::testDiscretization", "[PressureGradientPVE][TriP2P1P1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1()).testDiscretization();
}
TEST_CASE("PressureGradientPVE::TriP2P1P1::testResidual", "[PressureGradientPVE][TriP2P1P1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1()).testResidual();
}
TEST_CASE("PressureGradientPVE::TriP2P1P1::testJacobianTaylorSeries", "[PressureGradientPVE][TriP2P1P1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientPVE::TriP2P1P1::testJacobianFiniteDiff", "[PressureGradientPVE][TriP2P1P1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1()).testJacobianFiniteDiff();
}

// TriP3P2P2
TEST_CASE("PressureGradientPVE::TriP3P2P2::testDiscretization", "[PressureGradientPVE][TriP3P2P2][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2()).testDiscretization();
}
TEST_CASE("PressureGradientPVE::TriP3P2P2::testResidual", "[PressureGradientPVE][TriP3P2P2][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2()).testResidual();
}
TEST_CASE("PressureGradientPVE::TriP3P2P2::testJacobianTaylorSeries", "[PressureGradientPVE][TriP3P2P2][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientPVE::TriP3P2P2::testJacobianFiniteDiff", "[PressureGradientPVE][TriP3P2P2][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2()).testJacobianFiniteDiff();
}

// QuadQ2Q1Q1
TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1::testDiscretization", "[PressureGradientPVE][QuadQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1()).testDiscretization();
}
TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1::testResidual", "[PressureGradientPVE][QuadQ2Q1Q1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1()).testResidual();
}
TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1::testJacobianTaylorSeries", "[PressureGradientPVE][QuadQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1::testJacobianFiniteDiff", "[PressureGradientPVE][QuadQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1()).testJacobianFiniteDiff();
}

// QuadQ3Q2Q2
TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2::testDiscretization", "[PressureGradientPVE][QuadQ3Q2Q2][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2()).testDiscretization();
}
TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2::testResidual", "[PressureGradientPVE][QuadQ3Q2Q2][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2()).testResidual();
}
TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2::testJacobianTaylorSeries", "[PressureGradientPVE][QuadQ3Q2Q2][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2()).testJacobianTaylorSeries();
}
TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2::testJacobianFiniteDiff", "[PressureGradientPVE][QuadQ3Q2Q2][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2()).testJacobianFiniteDiff();
}

// TriP2P1P1 w/state variables
// TEST_CASE("PressureGradientPVE::TriP2P1P1_StateVars::testDiscretization", "[PressureGradientPVE][TriP2P1P1_StateVars][discretization]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1_StateVars()).testDiscretization();
// }
// TEST_CASE("PressureGradientPVE::TriP2P1P1_StateVars::testResidual", "[PressureGradientPVE][TriP2P1P1_StateVars][residual]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1_StateVars()).testResidual();
// }
// TEST_CASE("PressureGradientPVE::TriP2P1P1_StateVars::testJacobianTaylorSeries", "[PressureGradientPVE][TriP2P1P1_StateVars][Jacobian Taylor series]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1_StateVars()).testJacobianTaylorSeries();
// }
// TEST_CASE("PressureGradientPVE::TriP2P1P1_StateVars::testJacobianFiniteDiff", "[PressureGradientPVE][TriP2P1P1_StateVars][Jacobian finite difference]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP2P1P1_StateVars()).testJacobianFiniteDiff();
// }

// // TriP3P2P2 with state variables
// TEST_CASE("PressureGradientPVE::TriP3P2P2_StateVars::testDiscretization", "[PressureGradientPVE][TriP3P2P2_StateVars][discretization]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2_StateVars()).testDiscretization();
// }
// TEST_CASE("PressureGradientPVE::TriP3P2P2_StateVars::testResidual", "[PressureGradientPVE][TriP3P2P2_StateVars][residual]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2_StateVars()).testResidual();
// }
// TEST_CASE("PressureGradientPVE::TriP3P2P2_StateVars::testJacobianTaylorSeries", "[PressureGradientPVE][TriP3P2P2_StateVars][Jacobian Taylor series]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2_StateVars()).testJacobianTaylorSeries();
// }
// TEST_CASE("PressureGradientPVE::TriP3P2P2_StateVars::testJacobianFiniteDiff", "[PressureGradientPVE][TriP3P2P2_StateVars][Jacobian finite difference]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::TriP3P2P2_StateVars()).testJacobianFiniteDiff();
// }

// // QuadQ2Q1Q1 with state variables
// TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1_StateVars::testDiscretization", "[PressureGradientPVE][QuadQ2Q1Q1_StateVars][discretization]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1_StateVars()).testDiscretization();
// }
// TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1_StateVars::testResidual", "[PressureGradientPVE][QuadQ2Q1Q1_StateVars][residual]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1_StateVars()).testResidual();
// }
// TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1_StateVars::testJacobianTaylorSeries", "[PressureGradientPVE][QuadQ2Q1Q1_StateVars][Jacobian Taylor series]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1_StateVars()).testJacobianTaylorSeries();
// }
// TEST_CASE("PressureGradientPVE::QuadQ2Q1Q1_StateVars::testJacobianFiniteDiff", "[PressureGradientPVE][QuadQ2Q1Q1_StateVars][Jacobian finite difference]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ2Q1Q1_StateVars()).testJacobianFiniteDiff();
// }

// // QuadQ3Q2Q2 with state variables
// TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2_StateVars::testDiscretization", "[PressureGradientPVE][QuadQ3Q2Q2][discretization]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2_StateVars()).testDiscretization();
// }
// TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2_StateVars::testResidual", "[PressureGradientPVE][QuadQ3Q2Q2][residual]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2_StateVars()).testResidual();
// }
// TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2_StateVars::testJacobianTaylorSeries", "[PressureGradientPVE][QuadQ3Q2Q2][Jacobian Taylor series]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2_StateVars()).testJacobianTaylorSeries();
// }
// TEST_CASE("PressureGradientPVE::QuadQ3Q2Q2_StateVars::testJacobianFiniteDiff", "[PressureGradientPVE][QuadQ3Q2Q2][Jacobian finite difference]") {
//     pylith::TestLinearPoroviscoelasticity(pylith::PressureGradientPVE::QuadQ3Q2Q2_StateVars()).testJacobianFiniteDiff();
// }

// Poroviscoelasticity Terzaghi Problem
TEST_CASE("TerzaghiPVE::QuadQ2Q1Q1::testDiscretization", "[TerzaghiPVE][QuadQ2Q1Q1][discretization]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE::QuadQ2Q1Q1()).testDiscretization();
}
TEST_CASE("TerzaghiPVE::QuadQ2Q1Q1::testResidual", "[TerzaghiPVE][QuadQ2Q1Q1][residual]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE::QuadQ2Q1Q1()).testResidual();
}
TEST_CASE("TerzaghiPVE::QuadQ2Q1Q1::testJacobianTaylorSeries", "[TerzaghiPVE][QuadQ2Q1Q1][Jacobian Taylor series]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE::QuadQ2Q1Q1()).testJacobianTaylorSeries();
}
TEST_CASE("TerzaghiPVE::QuadQ2Q1Q1::testJacobianFiniteDiff", "[TerzaghiPVE][QuadQ2Q1Q1][Jacobian finite difference]") {
    pylith::TestLinearPoroviscoelasticity(pylith::TerzaghiPVE::QuadQ2Q1Q1()).testJacobianFiniteDiff();
}

// End of file
