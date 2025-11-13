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

#include "pylith/problems/InitialConditionDomain.hh" // implementation of class methods

#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/Field.hh" // USES Field

#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Constructor
pylith::problems::InitialConditionDomain::InitialConditionDomain(void) {
    PyreComponent::setName("initialconditiondomain");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::InitialConditionDomain::~InitialConditionDomain(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::InitialConditionDomain::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    _db = nullptr; // :KLLUDGE: Should use shared pointer.

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set spatial database holding initial conditions.
void
pylith::problems::InitialConditionDomain::setDB(const std::shared_ptr<spatialdata::spatialdb::SpatialDB>& db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setDB(db="<<db<<")");

    _db = db;

    PYLITH_METHOD_END;
} // setDB


// ---------------------------------------------------------------------------------------------------------------------
// Set solver type.
void
pylith::problems::InitialConditionDomain::setValues(const std::shared_ptr<pylith::topology::Field>& solution,
                                                    const spatialdata::units::Nondimensional& normalizer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setValues(solution="<<solution<<", normalizer)");
    assert(solution);

    pylith::topology::FieldQuery fieldQuery(solution);

    for (auto subfield : _subfields) {
        fieldQuery.addSubfield(subfield);
    } // for

    fieldQuery.open(_db, normalizer.getLengthScale());
    fieldQuery.query();
    fieldQuery.close();

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        PYLITH_COMPONENT_DEBUG("Displaying solution field");
        solution->view("Solution field with initial values", pylith::topology::Field::VIEW_ALL);
    } // if

    PYLITH_METHOD_END;
} // setValues


// End of file
