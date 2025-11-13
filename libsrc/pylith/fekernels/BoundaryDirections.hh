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
#include "pylith/utils/journals.hh"

#include <cassert> // USES assert()

class pylith::fekernels::BoundaryDirections {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /** Compute tangential directions for boundary in 3D.
     *
     * @param[out] tanDir1 First tangential direction.
     * @param[out] tanDIr2 Second tangential direction.
     * @param[in] refDir1 First choice for reference direction.
     * @param[in] refDir2 Second choice for reference direction if first fails.
     * @param[in] normDir Normal direction.
     */
    static inline
    void tangential_directions(pylith::real tanDir1[],
                               pylith::real tanDir2[],
                               const pylith::real refDir1[],
                               const pylith::real refDir2[],
                               const pylith::real normalDir[]) {
        assert(tanDir1);
        assert(tanDir2);
        assert(refDir1);
        assert(refDir2);
        assert(normalDir);

        const int dim = 3;

        // Choose reference direction 1 unless it nearly coincides with normal direction.
        pylith::real refDir[dim] = { refDir1[0], refDir1[1], refDir1[2] };
        if (fabs(refDir[0]*normalDir[0] + refDir[1]*normalDir[1] + refDir[2]*normalDir[2]) > 0.999) {
            for (int i = 0; i < dim; ++i) {
                refDir[i] = refDir2[i];
            } // for
        } // if

        // tanDir1 = refDir x normalDir
        _unitCrossProduct(tanDir1, refDir, normalDir);

        // tranDir2 = normalDir x tanDir1
        _unitCrossProduct(tanDir2, normalDir, tanDir1);
    }

    // ------------------------------------------------------------------------------------------------
    /** Normal direction.
     */
    static inline
    void normalDir(const pylith::integer dim,
                   const pylith::integer numS,
                   const pylith::integer numA,
                   const pylith::integer sOff[],
                   const pylith::integer sOff_x[],
                   const pylith::real s[],
                   const pylith::real s_t[],
                   const pylith::real s_x[],
                   const pylith::integer aOff[],
                   const pylith::integer aOff_x[],
                   const pylith::real a[],
                   const pylith::real a_t[],
                   const pylith::real a_x[],
                   const pylith::real t,
                   const pylith::real x[],
                   const pylith::real n[],
                   const pylith::integer numConstants,
                   const pylith::real constants[],
                   pylith::real normalDir[]) {
        const pylith::integer spaceDim = dim; // dim passed in is embedded dimension (spaceDim).

        for (pylith::integer i = 0; i < spaceDim; ++i) {
            normalDir[i] = n[i];
        } // for
    } // normalDir

    // ------------------------------------------------------------------------------------------------
    /** Horizontal tangential direction.
     *
     * If normal direction is parallel to the z-direction, then .
     */
    static inline
    void tangentialDirHoriz(const pylith::integer dim,
                            const pylith::integer numS,
                            const pylith::integer numA,
                            const pylith::integer sOff[],
                            const pylith::integer sOff_x[],
                            const pylith::real s[],
                            const pylith::real s_t[],
                            const pylith::real s_x[],
                            const pylith::integer aOff[],
                            const pylith::integer aOff_x[],
                            const pylith::real a[],
                            const pylith::real a_t[],
                            const pylith::real a_x[],
                            const pylith::real t,
                            const pylith::real x[],
                            const pylith::real n[],
                            const pylith::integer numConstants,
                            const pylith::real constants[],
                            pylith::real tangentialDir[]) {
        const pylith::integer spaceDim = dim; // dim passed in is embedded dimension (spaceDim).

        switch (spaceDim) {
        case 2: {
            const pylith::integer _dim = 2;
            const pylith::real tanDir[2] = {-n[1], n[0] };
            for (pylith::integer i = 0; i < _dim; ++i) {
                tangentialDir[i] = tanDir[i];
            } // for
            break;
        } // case 2
        case 3: {
            assert(numConstants >= 6);
            const pylith::integer _dim = 3;
            const pylith::real* refDir1 = &constants[0];
            const pylith::real* refDir2 = &constants[3];
            pylith::real tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

            for (pylith::integer i = 0; i < _dim; ++i) {
                tangentialDir[i] = tanDir1[i];
            } // for
            break;
        } // case 3
        default:
            assert(0);
        } // switch

    } // strikeDir

    // ------------------------------------------------------------------------------------------------
    /** Direction up dip (only valid in 3D).
     */
    static inline
    void tangentialDirVert(const pylith::integer dim,
                           const pylith::integer numS,
                           const pylith::integer numA,
                           const pylith::integer sOff[],
                           const pylith::integer sOff_x[],
                           const pylith::real s[],
                           const pylith::real s_t[],
                           const pylith::real s_x[],
                           const pylith::integer aOff[],
                           const pylith::integer aOff_x[],
                           const pylith::real a[],
                           const pylith::real a_t[],
                           const pylith::real a_x[],
                           const pylith::real t,
                           const pylith::real x[],
                           const pylith::real n[],
                           const pylith::integer numConstants,
                           const pylith::real constants[],
                           pylith::real tangentialDir[]) {
        const pylith::integer spaceDim = dim; // dim passed in is embedded dimension (spaceDim).

        switch (spaceDim) {
        case 2: {
            PYLITH_JOURNAL_LOGICERROR("Dip direction is not defined in 2D.");
            break;
        } // case 2
        case 3: {
            assert(numConstants >= 6);
            const pylith::integer _dim = 3;
            const pylith::real* refDir1 = &constants[0];
            const pylith::real* refDir2 = &constants[3];
            pylith::real tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

            for (pylith::integer i = 0; i < _dim; ++i) {
                tangentialDir[i] = tanDir2[i];
            } // for
            break;
        } // case 3
        default:
            assert(0);
        } // switch

    } // strikeDir

    /** Transform values from (tangential, normal) to (x, y).
     *
     * @param[out] valueXY Values in xy coordinate system.
     * @param[in] valuesTN Values in tagential-normal coordinate system.
     * @param[in] normalDir Normal direction unit vector.
     */
    static inline
    void toXY(pylith::real valuesXY[],
              const pylith::real valuesTN[],
              const pylith::real normalDir[]) {
        const pylith::integer _dim = 2;
        const pylith::real tanDir[2] = {-normalDir[1], normalDir[0] };

        for (pylith::integer i = 0; i < _dim; ++i) {
            valuesXY[i] += valuesTN[0]*tanDir[i] + valuesTN[1]*normalDir[i];
        } // for
    } // toXY

    /** Transform values from (tangential 1, tangential 2, normal) to (x, y, z).
     *
     * @param[out] valueXYZ Values in xyz coordinate system.
     * @param[in] valuesTN Values in tagential-normal coordinate system.
     * @param[in] refDir1 Reference direction 1.
     * @param[in] refDir2 Reference direction 2.
     * @param[in] normalDir Normal direction unit vector.
     */
    static inline
    void toXYZ(pylith::real valuesXYZ[],
               const pylith::real valuesTN[],
               const pylith::real refDir1[],
               const pylith::real refDir2[],
               const pylith::real normalDir[]) {
        const pylith::integer _dim = 3;

        pylith::real tanDir1[3], tanDir2[3];
        BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, normalDir);

        for (pylith::integer i = 0; i < _dim; ++i) {
            valuesXYZ[i] += valuesTN[0]*tanDir1[i] + valuesTN[1]*tanDir2[i] + valuesTN[2]*normalDir[i];
        } // for
    } // toXYZ

private:

    static inline
    void _unitCrossProduct(pylith::real result[],
                           const pylith::real vec1[],
                           const pylith::real vec2[]) {
        assert(result);
        assert(vec1);
        assert(vec2);

        result[0] = +vec1[1]*vec2[2] - vec1[2]*vec2[1];
        result[1] = +vec1[2]*vec2[0] - vec1[0]*vec2[2];
        result[2] = +vec1[0]*vec2[1] - vec1[1]*vec2[0];

        const int dim = 3;
        pylith::real resultMag = 0.0;
        for (int i = 0; i < dim; ++i) {
            resultMag += result[i]*result[i];
        } // for
        resultMag = sqrt(resultMag);
        assert(resultMag >= 1.0e-4);
        for (int i = 0; i < dim; ++i) {
            result[i] /= resultMag;
        } // for
    } // _unitCrossProduct

}; // BoundaryDirections

/* End of file */
