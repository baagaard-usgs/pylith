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

#include "pylith/topology/FieldBase.hh" // HASA validatorfn_type
#include "pylith/testing/testingfwd.hh" // USES FieldTester
#include "pylith/utils/utilsfwd.hh" // USES EventLogger

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys

#include <map> // HOLDSA std::map
#include <string> // HASA std::string

class pylith::topology::FieldQuery {
    friend class _FieldQuery;
    friend class TestFieldQuery; // unit testing
    friend class pylith::testing::FieldTester; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const pylith::real SCALE_TOLERANCE;

    // PUBLIC TYPEDEF /////////////////////////////////////////////////////////////////////////////
public:

    /** Function prototype for converter functions.
     *
     * @param[out] values Values for subfield.
     * @param[in] nvalues Number of values for subfield.
     * @param[in] Array of values from spatial database query.
     * @param[in] Indices of values from spatial database to use for computing subfield values.
     */
    typedef std::string (*convertfn_type)(pylith::scalar[],
                                          const pylith::integer,
                                          const pylith::scalar_array,
                                          const pylith::integer_array);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param field Field associated with query.
     */
    FieldQuery(const std::shared_ptr<pylith::topology::Field>& field);

    /// Destructor.
    ~FieldQuery(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Add subfield to query.
     *
     * Use the database passed in the call to openDB().
     *
     * @param[in] subfieldName Name of subfield.
     */
    void addSubfield(const std::string& subfieldName);

    /** Add subfield to query.
     *
     * Use the provided spatial database in the query.
     *
     * If the names of the values to query in the spatial database are not given via queryValues,
     * then the names of the subfield components are used.
     *
     * @param[in] subfieldName Name of subfield.
     * @param[in] db Spatial database to query (optional).
     * @param[in] queryValues Array of names of spatial database values for subfield.
     * @param[in] converter Function to convert spatial database values to subfield value (optional).
     */
    void addSubfield(const std::string& subfieldName,
                     std::shared_ptr<spatialdata::spatialdb::SpatialDB>& db,
                     const pylith::string_vector& queryValues=pylith::string_vector(),
                     convertfn_type converter=nullptr);

    /// Initialize query with default query information.
    void initializeWithDefaultQueries(void);

    /** Open spatial database query for setting values in field.
     *
     * @param db Spatial database to query.
     * @param lengthScale Length scale for dimensionalization of
     * location coordinates.
     */
    void open(const std::shared_ptr<spatialdata::spatialdb::SpatialDB>& db,
              const pylith::real lengthScale);

    /// Query spatial database to set values in field.
    void query(void);

    /** Query spatial database for points in label to set values in field.
     *
     * @param[in] labelName Name of label.
     * @param[in] labelValue Value of label.
     */
    void queryUsingLabel(const char* labelName,
                         const pylith::integer labelValue=1);

    /// Close spatial database query for setting values in field.
    void close(void);

    /** Query of values from spatial database point.
     *
     * Includes nondimensionalization but no conversion of values.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] t Current time.
     * @param[in] x Coordinates (nondimensioned) of point location for query.
     * @param[in] nvalues Size of values array.
     * @param[out] values Array of values to be returned.
     * @param[in] context Query context.
     * @returns PETSc error code (0 for success).
     */
    static
    PetscErrorCode queryPointFn(pylith::integer dim,
                                pylith::real t,
                                const pylith::real x[],
                                pylith::integer nvalues,
                                pylith::scalar* values,
                                void* context);

    // PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private:

    struct SubfieldQuery {
        std::shared_ptr<spatialdata::spatialdb::SpatialDB> db; ///< Spatial database to query.
        pylith::string_vector queryValues; ///< Values to use from spatial database.
        convertfn_type converter; ///< Function to convert spatial database values to subfield values.

        SubfieldQuery(void) :
            converter(nullptr) {}


    }; // SubfieldQuery

    typedef std::map<std::string, SubfieldQuery> subfieldquery_t;

    /// Function prototype for query functions.
    typedef PetscErrorCode (*queryfn_type)(pylith::integer,
                                           pylith::real,
                                           const pylith::real[],
                                           pylith::integer,
                                           pylith::scalar*,
                                           void*);

    /// Context for spatial database queries.
    struct DBQueryContext {
        std::shared_ptr<spatialdata::spatialdb::SpatialDB> db; ///< Spatial database.
        std::shared_ptr<spatialdata::geocoords::CoordSys> cs; ///< Coordinate system of input point locations.
        pylith::real lengthScale; ///< Length scale for dimensionalizing coordinates.
        pylith::real valueScale; ///< Scale for dimensionalizing values for subfield.
        pylith::real validatorTolerance; ///< Tolerance relative to valueScale for validation.
        std::string description; ///< Name of value;
        pylith::scalar_array queryValues; ///< Values returned by spatial database query;
        pylith::integer_array queryIndices; ///< Indices of spatial database values to use for subfield.
        convertfn_type converter; ///< Function to convert values to subfield (optional).
        pylith::topology::FieldBase::validatorfn_type validator; ///< Function to validate values (optional).

        DBQueryContext(void) :
            lengthScale(1.0),
            valueScale(1.0),
            validatorTolerance(0.0),
            description("unknown"),
            converter(nullptr),
            validator(nullptr) {}


    }; // DBQueryStruct

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    const std::shared_ptr<pylith::topology::Field>& _field; ///< Field associated with query.

    std::vector<queryfn_type> _functions; ///< Functions implementing queries.
    std::vector<DBQueryContext> _contexts; ///< Contexts for performing query for each subfield.
    std::vector<DBQueryContext*> _contextPtrs; ///< Contexts for performing query for each subfield.
    subfieldquery_t _subfieldQueries;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FieldQuery(const FieldQuery&) = delete;
    const FieldQuery& operator=(const FieldQuery&) = delete;

    // Validators /////////////////////////////////////////////////////////////////////////////////
public:

    class Validators {
public:

        /** Validator for positive values.
         *
         * If scale and tolerance are greater than zero, then the value must be
         * in the range [scale/tolerance, scale*tolerance].
         *
         * @param[in] value Value to validate.
         * @param[in] scale Scale for nondimensionalization.
         * @param[in] tolerance Tolerance relative to scale for validation.
         * @returns Error message if not positive, empty string otherwise.
         */
        static
        std::string positive(const pylith::real value,
                             const pylith::real scale,
                             const pylith::real tolerance);

        /** Validator for nonnegative values.
         *
         * If scale and tolerance are greater than zero, then the value must be
         * in the range [scale/tolerance, scale*tolerance].
         *
         * @param[in] value Value to validate.
         * @param[in] scale Scale for nondimensionalization.
         * @param[in] tolerance Tolerance relative to scale for validation.
         * @returns Error message if negative, empty string otherwise.
         */
        static
        std::string nonnegative(const pylith::real value,
                                const pylith::real scale,
                                const pylith::real tolerance);

        /** Validator for scale only.
         *
         * If scale and tolerance are greater than zero, then the value must be
         * in the range [scale/tolerance, scale*tolerance].
         *
         * @param[in] value Value to validate.
         * @param[in] scale Scale for nondimensionalization.
         * @param[in] tolerance Tolerance relative to scale for validation.
         * @returns Error message if negative, empty string otherwise.
         */
        static
        std::string scale(const pylith::real value,
                          const pylith::real scale,
                          const pylith::real tolerance);

    }; // Validators

}; // FieldQuery

// End of file
