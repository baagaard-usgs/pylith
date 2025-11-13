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

#include "pylith/problems/InitialCondition.hh" // implementation of class methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_

// ----------------------------------------------------------------------
// Constructor
pylith::problems::InitialCondition::InitialCondition(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::problems::InitialCondition::~InitialCondition(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::problems::InitialCondition::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Set fields for initial condition.
void
pylith::problems::InitialCondition::setSubfields(const pylith::string_vector& subfields) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("setFields(subfields size="<<subfields.size()<<")");

    _subfields = subfields;

    PYLITH_METHOD_END;
} // setFields


// ---------------------------------------------------------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::problems::InitialCondition::verifyConfiguration(const pylith::topology::Field& solution) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("verifyConfiguration(solution="<<solution.getName()<<")");

    for (auto subfield : _subfields) {
        if (!solution.hasSubfield(subfield.c_str())) {
            std::ostringstream msg;
            msg << "Cannot specify initial conditions for solution subfield '"<< subfield
                << "' in component '" << PyreComponent::getIdentifier() << "'"
                << "; field is not in solution.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    PYLITH_METHOD_END;
} // verityConfiguration


// End of file
