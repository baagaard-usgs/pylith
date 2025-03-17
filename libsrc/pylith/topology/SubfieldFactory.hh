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

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::Discretization
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery::convertfn_type

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Normalizer

class pylith::topology::SubfieldFactory : public pylith::utils::GenericComponent {
    friend class TestSubfieldFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SubfieldFactory(void);

    /// Destructor.
    virtual ~SubfieldFactory(void);

    /** Get number of subfield discretizations.
     *
     * @returns Number of subfield discretizations.
     */
    int getNumSubfields(void) const;

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] subfieldName Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] dimension Dimension of points for discretization.
     * @param[in] isFaultOnly True if subfield is limited to fault degrees of freedom.
     * @param[in] cellBasis Type of basis functions to use (e.g., simplex, tensor, or default).
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space.
     */
    void setSubfieldDiscretization(const char* subfieldName,
                                   const int basisOrder,
                                   const int quadOrder,
                                   const int dimension,
                                   const bool isFaultOnly,
                                   const pylith::topology::FieldBase::CellBasis cellBasis,
                                   const pylith::topology::FieldBase::SpaceEnum feSpace,
                                   const bool isBasisContinuous);

    /** Get discretization information for subfield.
     *
     * @param[in] subfieldName Name of subfield.
     * @return Discretization information for auxiliary subfield. If
     * discretization information was not set, then use "default".
     */
    const pylith::topology::FieldBase::Discretization& getSubfieldDiscretization(const char* subfieldName) const;

    /** Set spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void setQueryDB(std::shared_ptr<spatialdata::spatialdb::SpatialDB>& value);

    /** Get spatial database for filling auxiliary subfields.
     *
     * @returns Pointer to database.
     */
    const spatialdata::spatialdb::SpatialDB* getQueryDB(void) const;

    /** Initialize factory for setting up auxiliary subfields.
     *
     * @param[inout] field Field for which subfields are to be created.
     * @param[in] normalizer Scales for nondimensionalization.
     * @param[in] spaceDim Spatial dimension of problem.
     * @param[in] defaultDescription Default description for new subfields.
     */
    virtual
    void initialize(std::shared_ptr<pylith::topology::Field>& field,
                    const spatialdata::units::Nondimensional& normalizer,
                    const int spaceDim,
                    std::unique_ptr<pylith::topology::FieldBase::Description>& defaultDescription);

    /// Set subfield values using spatial database.
    void setValuesFromDB(void);

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

    std::shared_ptr<pylith::topology::Field> _field; ///< Field.
    pylith::topology::FieldBase::discretizations_map _subfieldDiscretizations; ///< Discretization for each subfield.

    /// Description for default subfield.
    std::unique_ptr<pylith::topology::FieldBase::Description> _defaultDescription;

    /// Database of values for filling subfields.
    std::shared_ptr<spatialdata::spatialdb::SpatialDB> _queryDB;

    /// Field query for filling subfield values via spatial database.
    std::unique_ptr<pylith::topology::FieldQuery> _fieldQuery;

    std::shared_ptr<spatialdata::units::Nondimensional> _normalizer; ///< Scales for nondimensionalization.
    size_t _spaceDim; ///< Spatial dimension.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    SubfieldFactory(const SubfieldFactory &); ///< Not implemented.
    const SubfieldFactory& operator=(const SubfieldFactory&); ///< Not implemented

}; // class SubfieldFactory

// End of file
