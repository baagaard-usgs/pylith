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
#include "MomentumLayout.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::momentum {
    template<int pde_dim, class StrainModel, class StressModel, class SolutionLayout, class AuxiliaryLayout>
    struct BodyForce;
} // namespce


/// f0 = f_i + ρ g_i
template<size_t dim, class StrainModel, class StressModel, class SolutionLayout, class AuxiliaryLayout>
struct pylith::fekernels::momentum::BodyForce {
    PYLITH_KERNEL static void f0(const pylith::integer cellDim,
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
                                 pylith::scalar f0[]) noexcept {
        assert(dim == cellDim);

        const auto auxiliary = AuxiliaryLayout::template unpack<dim>(aOff, aOff_x, a, a_t, a_x);

        if constexpr (AuxiliaryLayout::has(MOMENTUM_BODY_FORCE)) {
            const auto& body_force = auxiliary.template get<MOMENTUM_BODY_FORCE>();
            for (size_t i = 0; i < dim; i++) {
                f0[i] += body_force(i);
            } // for
        } // if

        if constexpr (AuxiliaryLayout::has(MOMENTUM_GRAVITY)) {
            const auto& density = auxiliary.density;
            const auto& grav_acc = auxiliary.template get<MOMENTUM_GRAVITY>();
            for (size_t i = 0; i < dim; i++) {
                f0[i] += density * grav_acc(i);
            } // for
        } // if

    } // f0

}; // BodyForce
