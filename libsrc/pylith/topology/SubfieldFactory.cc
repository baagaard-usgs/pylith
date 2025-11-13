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

#include "pylith/topology/SubfieldFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::topology::SubfieldFactory::SubfieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::SubfieldFactory::~SubfieldFactory(void) {}


// ------------------------------------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::topology::SubfieldFactory::setDiscretization(const char* subfieldName,
                                                     const pylith::topology::FieldBase::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setDiscretization(subfieldName="<<subfieldName<<")");
    assert(discretization.dimension != 0);

    _subfields[subfieldName] = discretization;

    PYLITH_METHOD_END;
} // setDiscretization


// ------------------------------------------------------------------------------------------------
// Check if factory has discretization information for subfield.
bool
pylith::topology::SubfieldFactory::hasDiscretization(const char* subfieldName) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfield="<<subfieldName<<")");

    subfields_t::const_iterator iter = _subfields.find(subfieldName);
    bool found = iter != _subfields.end();

    PYLITH_METHOD_RETURN(found);
} // hasDiscretization


// ------------------------------------------------------------------------------------------------
// Add subfield to field.
void
pylith::topology::SubfieldFactory::addSubfield(const std::string& subfieldName) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addSubfield="<<subfieldName<<")");

    subfields_t::const_iterator iter = _subfields.find(subfieldName);
    if (iter == _subfields.end()) {
        iter = _subfields.find("default");
    } // if
    if (iter == _subfields.end()) {
        PYLITH_JOURNAL_LOGICERROR("Could not find subfield " << subfieldName << " in subfield factory.");
    } // if

    const size_t spaceDim = _field->getSpaceDim();
    const pylith::topology::FieldBase::Description& description = _getDescription(subfieldName, spaceDim);
    _field->subfieldAdd(description, iter->second);

    PYLITH_METHOD_END;
} // addSubfield


// ------------------------------------------------------------------------------------------------
// Initialize factory for setting up auxiliary subfields.
void
pylith::topology::SubfieldFactory::open(std::shared_ptr<pylith::topology::Field>& field,
                                        const std::shared_ptr<spatialdata::units::Nondimensional>& normalizer) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("open(field="<<field<<", normalizer="<<&normalizer<<")");
    assert(field);

    _field = field;
    _normalizer = normalizer;
    assert(2 == field->getSpaceDim() || 3 == field->getSpaceDim());

    PYLITH_METHOD_END;
} // open


// ------------------------------------------------------------------------------------------------
// Close factory after setting up auxiliary subfields.
void
pylith::topology::SubfieldFactory::close(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("close()");

    _field.reset();
    _normalizer.reset();

    PYLITH_METHOD_END;
} // close


// End of file
