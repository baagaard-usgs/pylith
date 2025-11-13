/*
 * ================================================================================================
 * This code is part of PyLith, developed through the Computational Infrastructure
 * for Geodynamics (https://github.com/geodynamics/pylith).
 *
 * Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
 * All rights reserved.
 *
 * See https://mit-license.org/ and LICENSE.md and for license information.
 * =================================================================================================
 */
#pragma once

#include "pylith/fekernels/fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::Solution {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Kernel interface.
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
     * @param[out] f0 [dim].
     */

    /** Pass thru solution in pointwise function.
     *
     * We pass through the solution to the resulting field. The auxiliary field is ignored.
     */
    static inline
    void passThruSubfield(const pylith::integer dim,
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
                          pylith::scalar field[]) {
        assert(s);
        assert(sOff);
        assert(field);
        const pylith::integer subfieldIndex = pylith::integer(t); // :KLUDGE: Easiest way to get subfield to extract
                                                                  // into fn.

        const pylith::integer sEnd = sOff[subfieldIndex+1];
        for (pylith::integer iS = sOff[subfieldIndex], iF = 0; iS < sEnd; ++iS, ++iF) {
            field[iF] = s[iS];
        } // for
    } // passThruSubfield

}; // Solution

// End of file
