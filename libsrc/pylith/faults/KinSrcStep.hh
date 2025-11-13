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

/** KinSrcStep
 *
 * @brief Step slip-time function.
 *
 * Slip time function is a step function at the time of slip.
 *
 * slip = final_slip for t >= t0.
 *
 * slip = 0 for t < t0.
 */
class pylith::faults::KinSrcStep : public KinSrc {
    friend class TestKinSrcStep; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcStep(void);

    /// Destructor.
    ~KinSrcStep(void);

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

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                              const spatialdata::geocoords::CoordSys* cs);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    KinSrcStep(const KinSrcStep&) = delete;
    const KinSrcStep& operator=(const KinSrcStep&) = delete;

}; // class KinSrcStep

// End of file
