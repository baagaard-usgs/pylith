// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations
#include "pylith/topology/FieldFactory.hh" // ISA FieldFactory

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery::queryfn_type

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Normalizer

class pylith::feassemble::AuxiliaryFactory : public pylith::topology::FieldFactory {
    friend class TestAuxiliaryFactory; // unit testing

    // PUBLIC CONTS ///////////////////////////////////////////////////////////////////////////////
public:

    static const PylithReal SCALE_TOLERANCE; ///< Tolerance for validating nondimensionalized values.

    // PUBLIC TYPEDEGS ////////////////////////////////////////////////////////////////////////////
public:

    typedef void (*validatorfn_type)(const pylith::topology::Field* /* field */);

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactory(void);

    /// Destructor.
    virtual ~AuxiliaryFactory(void);

    /** Set spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void setQueryDB(spatialdata::spatialdb::SpatialDB* value);

    /** Get spatial database for filling auxiliary subfields.
     *
     * @returns Pointer to database.
     */
    const spatialdata::spatialdb::SpatialDB* getQueryDB(void) const;

    /** Initialize factory for setting up auxiliary subfields.
     *
     * @param[inout] field Auxiliary field for which subfields are to be created.
     * @param[in] normalizer Scales for nondimensionalization.
     * @param[in] spaceDim Spatial dimension of problem.
     * @param[in] defaultDescription Default description for new subfields.
     */
    void initialize(pylith::topology::Field* field,
                    const spatialdata::units::Nondimensional& normalizer,
                    const int spaceDim,
                    const pylith::topology::FieldBase::Description* defaultDescription=NULL);

    /** Set subfield values using spatial database.
     *
     * @param[in] validateFn Function for validating auxiliary field.
     */
    void setValuesFromDB(validatorfn_type valifateFn=NULL);

    /** Set query function for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @param[in] namesDBValues Array of names of values to use from spatial database.
     * @param[in] numDBValues Size of names array.
     * @param[in] convertFn Function to convert spatial database values to subfield values.
     * @param[in] db Spatial database to query.
     */
    void setSubfieldQuery(const char* subfieldName,
                          const char* namesDBValues[]=NULL,
                          const size_t numDBValues=0,
                          pylith::topology::FieldQuery::convertfn_type convertFn=NULL,
                          spatialdata::spatialdb::SpatialDB* db=NULL);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Database of values for filling subfields.
     * Auxiliary subfields are not filled if NULL.
     *
     * Currently, this is the only way to populate the auxiliary subfields.
     */
    spatialdata::spatialdb::SpatialDB* _queryDB;

    /// Field query for filling subfield values via spatial database.
    pylith::topology::FieldQuery* _fieldQuery;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactory(const AuxiliaryFactory &); ///< Not implemented.
    const AuxiliaryFactory& operator=(const AuxiliaryFactory&); ///< Not implemented

}; // class AuxiliaryFactory

// End of file
