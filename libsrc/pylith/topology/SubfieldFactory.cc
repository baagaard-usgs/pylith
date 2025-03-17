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

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/error.hh" // USES PYLITH_METHOD*
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <utility> // USES std::move
#include <cassert>

// ------------------------------------------------------------------------------------------------
// Default constructor.
pylith::topology::SubfieldFactory::SubfieldFactory(void) :
    _normalizer(new spatialdata::units::Nondimensional),
    _spaceDim(0) {
    GenericComponent::setName("auxiliaryfactory");
    _subfieldDiscretizations["default"] = pylith::topology::FieldBase::Discretization();
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::SubfieldFactory::~SubfieldFactory(void) {
    _field.reset();
    _defaultDescription.reset();
    _normalizer.reset();
    _queryDB.reset();
    _fieldQuery.reset();
} // destructor


// ------------------------------------------------------------------------------------------------
// Get number of subfield discretizations.
int
pylith::topology::SubfieldFactory::getNumSubfields(void) const {
    return _subfieldDiscretizations.size();
} // getNumSubfields


// ------------------------------------------------------------------------------------------------
// Set discretization information for auxiliary subfield.
void
pylith::topology::SubfieldFactory::setSubfieldDiscretization(const char* subfieldName,
                                                             const int basisOrder,
                                                             const int quadOrder,
                                                             const int dimension,
                                                             const bool isFaultOnly,
                                                             const pylith::topology::FieldBase::CellBasis cellBasis,
                                                             const pylith::topology::FieldBase::SpaceEnum feSpace,
                                                             const bool isBasisContinuous) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSubfieldDiscretization(subfieldName="<<subfieldName<<", basisOrder="<<basisOrder<<", quadOrder="<<quadOrder<<", dimension="<<dimension<<", cellBasis="<<cellBasis<<", isBasisContinuous="<<isBasisContinuous<<")");
    assert(dimension != 0);

    pylith::topology::FieldBase::Discretization discretization;
    discretization.basisOrder = basisOrder;
    discretization.quadOrder = quadOrder;
    discretization.dimension = dimension;
    discretization.cellBasis = cellBasis;
    discretization.isFaultOnly = isFaultOnly;
    discretization.feSpace = feSpace;
    discretization.isBasisContinuous = isBasisContinuous;
    _subfieldDiscretizations[subfieldName] = discretization;

    PYLITH_METHOD_END;
} // setSubfieldDiscretization


// ------------------------------------------------------------------------------------------------
// Get discretization information for subfield.
const pylith::topology::FieldBase::Discretization&
pylith::topology::SubfieldFactory::getSubfieldDiscretization(const char* subfieldName) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("getSubfieldDiscretization(subfieldName="<<subfieldName<<")");

    pylith::topology::FieldBase::discretizations_map::const_iterator iter = _subfieldDiscretizations.find(subfieldName);
    if (iter != _subfieldDiscretizations.end()) {
        PYLITH_METHOD_RETURN(iter->second);
    } else { // not found so try default
        iter = _subfieldDiscretizations.find("default");
        if (iter == _subfieldDiscretizations.end()) {
            throw std::logic_error("Default discretization not set in field factory.");
        } // if
    } // if/else

    PYLITH_METHOD_RETURN(iter->second); // default
} // getSubfieldDiscretization


// ------------------------------------------------------------------------------------------------
// Set database for filling auxiliary subfields.
void
pylith::topology::SubfieldFactory::setQueryDB(std::shared_ptr<spatialdata::spatialdb::SpatialDB>& value) {
    _queryDB = value;
} // setQueryDB


// ------------------------------------------------------------------------------------------------
// Get database for filling auxiliary subfields.
const spatialdata::spatialdb::SpatialDB*
pylith::topology::SubfieldFactory::getQueryDB(void) const {
    return _queryDB.get();
} // getQueryDB


// ------------------------------------------------------------------------------------------------
// Initialize factory for setting up auxiliary subfields.
void
pylith::topology::SubfieldFactory::initialize(std::shared_ptr<pylith::topology::Field>& field,
                                              const spatialdata::units::Nondimensional& normalizer,
                                              const int spaceDim,
                                              std::unique_ptr<pylith::topology::FieldBase::Description>& defaultDescription) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("initialize(field="<<field<<", normalizer="<<&normalizer<<", spaceDim="<<spaceDim<<", defaultDescription="<<defaultDescription.get()<<")");

    assert(field);

    _field = field;
    if (defaultDescription) {
        _defaultDescription = std::move(defaultDescription);
    } else {
        _defaultDescription.reset();
    } // if/else
    assert(_normalizer);
    *_normalizer = normalizer;
    _spaceDim = spaceDim;

    assert(1 <= _spaceDim && _spaceDim <= 3);
    std::unique_ptr<pylith::topology::FieldQuery> tmpQuery(new pylith::topology::FieldQuery(*_field.get()));
    _fieldQuery = std::move(tmpQuery);

    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Initialize subfields.
void
pylith::topology::SubfieldFactory::setValuesFromDB(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setValuesFromDB()");

    assert(_normalizer);
    assert(_field);

    if (_queryDB) {
        assert(_fieldQuery);
        _fieldQuery->openDB(_queryDB.get(), _normalizer->getLengthScale());
        _fieldQuery->queryDB();
        _fieldQuery->closeDB(_queryDB.get());
    } else {
        PYLITH_JOURNAL_ERROR("Unknown case for filling auxiliary subfields.");
        throw std::logic_error("Unknown case for filling auxiliary subfields.");
    } // if/else

    _fieldQuery.reset();
    _field.reset();

    PYLITH_METHOD_END;
} // setValuesFromDB


// ------------------------------------------------------------------------------------------------
// Set query function for subfield.
void
pylith::topology::SubfieldFactory::setSubfieldQuery(const char* subfieldName,
                                                    const char* namesDBValues[],
                                                    const size_t numDBValues,
                                                    pylith::topology::FieldQuery::convertfn_type convertFn,
                                                    spatialdata::spatialdb::SpatialDB* db) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("setSubfieldQuery(subfieldName="<<subfieldName<<", namesDBValues="<<namesDBValues<<", numDBValues="<<numDBValues<<", convertFn="<<convertFn<<", db="<<db<<")");

    assert(_fieldQuery);
    _fieldQuery->setQuery(subfieldName, namesDBValues, numDBValues, convertFn, db);

    PYLITH_METHOD_END;
} // _setSubfieldQueryFn


// End of file
