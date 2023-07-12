// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "catch2/catch_test_macros.hpp"

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace utils {
        class TestPetscStack;
    }
}

class pylith::utils::TestPetscStack : public pylith::utils::GenericComponent {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    TestPetscStack(const bool missingBegin0,
                   const bool missingReturn0,
                   const bool missingBegin1,
                   const bool missingReturn1);

    /// Test detection of mismatches in PetscFunctionBeginUser() and PetscFunctionReturn()
    static
    void testStack(void);

    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////
private:

    PetscErrorCode _layer0(void);

    PetscErrorCode _layer1(void);

    static
    void _checkError(PetscErrorCode err);

    // PRIVATE MEMEBRS ////////////////////////////////////////////////////////////////////////////
private:

    bool _hasBegin0;
    bool _hasReturn0;
    bool _hasBegin1;
    bool _hasReturn1;

}; // class TestPetscStack

// ------------------------------------------------------------------------------------------------
TEST_CASE("TestPetscStack::testStack", "[TestPetscStack]") {
    pylith::utils::TestPetscStack::testStack();
}

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::utils::TestPetscStack::TestPetscStack(const bool hasBegin0,
                                              const bool hasReturn0,
                                              const bool hasBegin1,
                                              const bool hasReturn1) :
    _hasBegin0(hasBegin0),
    _hasReturn0(hasReturn0),
    _hasBegin1(hasBegin1),
    _hasReturn1(hasReturn1) {}


// ------------------------------------------------------------------------------------------------
void
pylith::utils::TestPetscStack::testStack(void) {
    PetscErrorCode err = PETSC_SUCCESS;

    SECTION("Correct stack") {
        TestPetscStack stack(true, true, true, true);
        err = stack._layer0();REQUIRE_NOTHROW(_checkError(err));
    }

    SECTION("Missing begin0") {
        TestPetscStack stack(false, true, true, true);
        // err = stack._layer0();REQUIRE_THROWS_AS(_checkError(err), std::runtime_error);
    }

    SECTION("Missing end0") {
        TestPetscStack stack(true, false, true, true);
        err = stack._layer0();REQUIRE_THROWS_AS(_checkError(err), std::runtime_error);
    }

    SECTION("Missing begin1") {
        TestPetscStack stack(true, true, false, true);
        err = stack._layer0();REQUIRE_THROWS_AS(_checkError(err), std::runtime_error);
    }

    SECTION("Missing return1") {
        TestPetscStack stack(true, true, true, false);
        err = stack._layer0();REQUIRE_THROWS_AS(_checkError(err), std::runtime_error);
    }

    SECTION("Missing begin/return 0") {
        TestPetscStack stack(false, false, true, true);
        err = stack._layer0();REQUIRE_NOTHROW(_checkError(err));
    }

    SECTION("Missing begin/return 1") {
        TestPetscStack stack(true, true, false, false);
        err = stack._layer0();REQUIRE_NOTHROW(_checkError(err));
    }

    SECTION("Missing begin 0 and return 1") {
        TestPetscStack stack(false, true, true, false);
        err = stack._layer0();REQUIRE_THROWS_AS(_checkError(err), std::runtime_error);
    }

    SECTION("Missing return 0 and begin 1") {
        TestPetscStack stack(true, false, false, true);
        err = stack._layer0();REQUIRE_THROWS_AS(_checkError(err), std::runtime_error);
    }

}


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::utils::TestPetscStack::_layer0(void) {
    if (_hasBegin0) {
        PYLITH_METHOD_BEGIN;
    }

    PetscErrorCode err = PETSC_SUCCESS;
    err = _layer1();PYLITH_CHECK_ERROR(err);

    if (_hasReturn0) {
        PYLITH_METHOD_RETURN(err);
    } else {
        return err;
    }
}


// ------------------------------------------------------------------------------------------------
PetscErrorCode
pylith::utils::TestPetscStack::_layer1(void) {
    if (_hasBegin1) {
        PYLITH_METHOD_BEGIN;
    }

    PetscErrorCode err = PETSC_SUCCESS;

    if (_hasReturn1) {
        PYLITH_METHOD_RETURN(err);
    } else {
        return err;
    }
}


// ------------------------------------------------------------------------------------------------
void
pylith::utils::TestPetscStack::_checkError(PetscErrorCode err) {
    PYLITH_CHECK_ERROR(err);
}
