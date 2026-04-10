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

#include "pylith/fekernels/common/kernel.hh"
#include "pylith/fekernels/common/Fields.hh"
#include "pylith/fekernels/common/Matrix.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::momentum {
    template<size_t im, class StrainModel, class StressModel, class SolutionLayout, class AuxiliaryLayout>
    struct DivStress;
} // pylith::fekernels::momentum


template<size_t dim, class StrainModel, class StressModel, class SolutionLayout, class AuxiliaryLayout>
struct pylith::fekernels::momentum::DivStress {
    // f1 = -σ
    PYLITH_KERNEL static void f1(const pylith::integer cellDim,
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
                                 const pylith::real x[],
                                 const pylith::integer numConstants,
                                 const pylith::scalar constants[],
                                 pylith::scalar f1[]) noexcept {
        assert(dim == cellDim);

        pylith::fekernels::common::Matrix<dim> strain;
        pylith::fekernels::common::Matrix<dim> stress;

        const auto solution = SolutionLayout::template unpack<dim>(sOff, sOff_x, s, s_t, s_x);
        const auto auxiliary = AuxiliaryLayout::template unpack<dim>(aOff, aOff_x, a, a_t, a_x);

        StrainModel::compute(strain, solution);
        StressModel::cauchyStress(stress, strain, auxiliary);

        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                f1[i*dim+j] = -stress(i, j);
            } // for
        } // for
    } // f1

    // Jf3uu = j(f,g,df,dg) = C(f,df,g,dg)
    PYLITH_KERNEL static void Jf3uu(const pylith::integer cellDim,
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
                                    const pylith::real x[],
                                    const pylith::integer numConstants,
                                    const pylith::scalar constants[],
                                    pylith::scalar Jf3[]) noexcept {
        assert(dim == cellDim);

        const auto auxiliary = AuxiliaryLayout::template unpack<dim>(aOff, aOff_x, a, a_t, a_x);
        constexpr pylith::scalar sign = -1.0;
        StressModel::cauchyStressTangent(Jf3, auxiliary, sign);
    } // Jf3

}; // DivStress
