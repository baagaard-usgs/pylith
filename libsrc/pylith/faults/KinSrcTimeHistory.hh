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

#include "pylith/faults/KinSrc.hh"

class pylith::faults::KinSrcTimeHistory : public KinSrc {
    friend class TestKinSrcTimeHistory; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcTimeHistory(void);

    /// Destructor.
    ~KinSrcTimeHistory(void);

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory* getTimeHistoryDB(void);

    /** Get requested slip subfields at time t.
     *
     * @param[inout] slipLocalVec Local PETSc vector for slip, slip rate, or slip acceleration values.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     * @param[in] bitSlipSubfields Slip subfields to compute.
     */
    void getSlipSubfields(PetscVec slipLocalVec,
                          pylith::topology::Field* faultAuxiliaryField,
                          const pylith::real t,
                          const pylith::real timeScale,
                          const int bitSlipSubfields);

    /** Slip time function kernel.
     *
     * The "solution" field s is ignored.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] numS Number of registered subfields in solution field.
     * @param[in] numA Number of registered subfields in auxiliary field.
     * @param[in] sOff Offset of registered subfields in solution field [numS].
     * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
     * @param[in] s Solution field with all subfields.
     * @param[in] s_t Time derivative of solution field.
     * @param[in] s_x Gradient of solution field.
     * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
     * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
     * @param[in] a Auxiliary field with all subfields.
     * @param[in] a_t Time derivative of auxiliary field.
     * @param[in] a_x Gradient of auxiliary field.
     * @param[in] t Time for residual evaluation.
     * @param[in] x Coordinates of point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] slip [dim].
     */
    static
    void slipFn(const pylith::integer dim,
                const pylith::integer numS,
                const pylith::integer numA,
                const pylith::integer sOff[],
                const pylith::integer sOff_x[],
                const pylith::scalar s[],
                const pylith::scalar s_t[],
                const pylith::scalar s_x[],
                const pylith::integer aOff[],
                const pylith::integer aOff_x[],
                const pylith::scalar a[],
                const pylith::scalar a_t[],
                const pylith::scalar a_x[],
                const pylith::real t,
                const pylith::scalar x[],
                const pylith::integer numConstants,
                const pylith::scalar constants[],
                pylith::scalar slip[]);

    /** Slip rate time function kernel.
     *
     * The "solution" field s is ignored.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] numS Number of registered subfields in solution field.
     * @param[in] numA Number of registered subfields in auxiliary field.
     * @param[in] sOff Offset of registered subfields in solution field [numS].
     * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
     * @param[in] s Solution field with all subfields.
     * @param[in] s_t Time derivative of solution field.
     * @param[in] s_x Gradient of solution field.
     * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
     * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
     * @param[in] a Auxiliary field with all subfields.
     * @param[in] a_t Time derivative of auxiliary field.
     * @param[in] a_x Gradient of auxiliary field.
     * @param[in] t Time for residual evaluation.
     * @param[in] x Coordinates of point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] slipRate [dim].
     */
    static
    void slipRateFn(const pylith::integer dim,
                    const pylith::integer numS,
                    const pylith::integer numA,
                    const pylith::integer sOff[],
                    const pylith::integer sOff_x[],
                    const pylith::scalar s[],
                    const pylith::scalar s_t[],
                    const pylith::scalar s_x[],
                    const pylith::integer aOff[],
                    const pylith::integer aOff_x[],
                    const pylith::scalar a[],
                    const pylith::scalar a_t[],
                    const pylith::scalar a_x[],
                    const pylith::real t,
                    const pylith::scalar x[],
                    const pylith::integer numConstants,
                    const pylith::scalar constants[],
                    pylith::scalar slipRate[]);

    /** Slip acceleration time function kernel.
     *
     * The "solution" field s is ignored.
     *
     * @param[in] dim Spatial dimension.
     * @param[in] numS Number of registered subfields in solution field.
     * @param[in] numA Number of registered subfields in auxiliary field.
     * @param[in] sOff Offset of registered subfields in solution field [numS].
     * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
     * @param[in] s Solution field with all subfields.
     * @param[in] s_t Time derivative of solution field.
     * @param[in] s_x Gradient of solution field.
     * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
     * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
     * @param[in] a Auxiliary field with all subfields.
     * @param[in] a_t Time derivative of auxiliary field.
     * @param[in] a_x Gradient of auxiliary field.
     * @param[in] t Time for residual evaluation.
     * @param[in] x Coordinates of point evaluation.
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] slipAcc [dim].
     */
    static
    void slipAccFn(const pylith::integer dim,
                   const pylith::integer numS,
                   const pylith::integer numA,
                   const pylith::integer sOff[],
                   const pylith::integer sOff_x[],
                   const pylith::scalar s[],
                   const pylith::scalar s_t[],
                   const pylith::scalar s_x[],
                   const pylith::integer aOff[],
                   const pylith::integer aOff_x[],
                   const pylith::scalar a[],
                   const pylith::scalar a_t[],
                   const pylith::scalar a_x[],
                   const pylith::real t,
                   const pylith::scalar x[],
                   const pylith::integer numConstants,
                   const pylith::scalar constants[],
                   pylith::scalar slipAcc[]);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                              const spatialdata::geocoords::CoordSys* cs);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    spatialdata::spatialdb::TimeHistory* _dbTimeHistory; ///< Time history database.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    KinSrcTimeHistory(const KinSrcTimeHistory&) = delete;
    const KinSrcTimeHistory& operator=(const KinSrcTimeHistory&) = delete;

}; // class KinSrcTimeHistory

// End of file
