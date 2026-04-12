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

#include "pylith/utils/types.hh"

#include "pylith/fekernels/common/kernel.hh"
#include "pylith/fekernels/common/ArgFields.hh"
#include "pylith/fekernels/common/OptionalFields.hh"

#include <cassert>
#include <cstddef>


namespace pylith::fekernels::pde::elasticity {
    // Flags for elasticity solution
    enum SolutionFlags : size_t {
        DEFAULT=0,
        FAULT=1 << 0,
        INERTIA=1 << 1,
    }; // SolutionFlags

    template<SolutionFlags flags> struct SolutionLayout;
} // namespace

/// Layout of solution field for elasticity
template<pylith::fekernels::pde::elasticity::SolutionFlags flags>
struct pylith::fekernels::pde::elasticity::SolutionLayout {
    // Is a given flag present?
    static constexpr bool has(SolutionFlags f) {
        return (flags & f);
    } // has

    // Use shared OptionalMember from common utilities
    using OptionalMember = pylith::fekernels::common::OptionalMember;

    /// Order of solution subfields
    enum Fields : int {
        DISPLACEMENT=0,
        VELOCITY=1,
        LAGRANGE_MULTIPLIER_FAULT=0 + pylith::fekernels::common::addIfPresent(has(INERTIA)),
        NUM_FIELDS=1 + pylith::fekernels::common::addIfPresent(has(INERTIA))
                    + pylith::fekernels::common::addIfPresent(has(FAULT))
    };

    /// Struct wth names for holding subfields
    template <size_t dim>
    struct Unpacked {
        pylith::fekernels::common::VectorField<dim> displacement;

        // Optional — zero size if absent
        [[no_unique_address]]
        OptionalMember<has(INERTIA), pylith::fekernels::common::VectorField<dim> > velocity;

        [[no_unique_address]]
        OptionalMember<has(FAULT), pylith::fekernels::common::VectorField<dim> > lagrange_multiplier_fault;

        // Type-safe accessor — compile error if field not present
        template <SolutionFlags F>
        PYLITH_KERNEL auto& get() {
            static_assert(has(F), "Elasticity solution subfield not present in this layout.");
            if constexpr (F == INERTIA) { return velocity.member;}
            if constexpr (F == FAULT) { return lagrange_multiplier_fault.member;}
        } // get()

    }; // Unpacked

    /// Unpack solution fields from array into struct with names.
    template <size_t dim>
    PYLITH_KERNEL static Unpacked<dim> unpack(const pylith::integer sOff[],
                                              const pylith::integer sOff_x[],
                                              const pylith::scalar s[],
                                              const pylith::scalar s_t[],
                                              const pylith::scalar s_x[]) {
        Unpacked<dim> data;

        data.displacement = {
            &s[sOff[DISPLACEMENT]],
            s_t ? &s_t[sOff[DISPLACEMENT]] : nullptr,
            s_x ? &s_x[sOff_x[DISPLACEMENT]] : nullptr,
        };

        if constexpr (has(INERTIA)) {
            data.velocity = {
                &s[sOff[VELOCITY]],
                s_t ? &s_t[sOff[VELOCITY]] : nullptr,
                s_x ? &s_x[sOff_x[VELOCITY]] : nullptr,
            };
        }

        if constexpr (has(FAULT)) {
            data.lagrange_multiplier_fault = {
                &s[sOff[LAGRANGE_MULTIPLIER_FAULT]],
                s_t ? &s_t[sOff[LAGRANGE_MULTIPLIER_FAULT]] : nullptr,
                s_x ? &s_x[sOff_x[LAGRANGE_MULTIPLIER_FAULT]] : nullptr,
            };
        }

        return data;
    } // unpack

}; // SolutionLayout
