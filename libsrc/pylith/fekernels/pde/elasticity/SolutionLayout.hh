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
#include "pylith/fekernels/common/Utils.hh"
#include "pylith/fekernels/common/ArgFields.hh"
#include "pylith/fekernels/common/OptionalFields.hh"

#include <cassert>
#include <cstddef>


namespace common = pylith::fekernels::common;

namespace pylith::fekernels::pde::elasticity {
    // Flags for elasticity solution
    enum class SolutionFlags : size_t {
        DEFAULT=0,
        FAULT=1 << 0,
        INERTIA=1 << 1,
    }; // SolutionFlags

    PYLITH_KERNEL constexpr SolutionFlags
    operator|(SolutionFlags a,
              SolutionFlags b) noexcept {
        return static_cast<SolutionFlags>(
            static_cast<size_t>(a) | static_cast<size_t>(b));
    } // operator|


    PYLITH_KERNEL constexpr SolutionFlags&
    operator|=(SolutionFlags& a,
               SolutionFlags b) noexcept {
        a = a | b;
        return a;
    }


    PYLITH_KERNEL constexpr bool
    operator&(SolutionFlags a,
              SolutionFlags b) noexcept {
        return static_cast<size_t>(a) & static_cast<size_t>(b);
    }


    template<SolutionFlags flags> struct SolutionLayout;
} // namespace


/// Layout of solution field for elasticity
template<pylith::fekernels::pde::elasticity::SolutionFlags flags>
struct pylith::fekernels::pde::elasticity::SolutionLayout {
    // Is a given flag present?
    static constexpr bool has(SolutionFlags f) {
        return (flags & f);
    } // has

    /// Order of solution subfields
    enum class Fields : uint32_t {
        DISPLACEMENT=0,
        INERTIA=1,
        LAGRANGE_MULTIPLIER_FAULT=1 + pylith::fekernels::common::addIf(has(SolutionFlags::INERTIA)),
        NUM_FIELDS=1 + pylith::fekernels::common::addIf(has(SolutionFlags::INERTIA))
                    + pylith::fekernels::common::addIf(has(SolutionFlags::FAULT))
    };

    /// Struct wth names for holding subfields
    template <typename Dim>
    struct Unpacked {
        using Layout = pylith::fekernels::pde::elasticity::SolutionLayout<flags>;

        pylith::fekernels::common::VectorField<Dim> displacement;

        // Optional — zero size if absent
        [[no_unique_address]]
        pylith::fekernels::common::OptionalMember<has(SolutionFlags::INERTIA), pylith::fekernels::common::VectorField<Dim> > velocity;

        [[no_unique_address]]
        pylith::fekernels::common::OptionalMember<has(SolutionFlags::FAULT), pylith::fekernels::common::VectorField<Dim> > lagrange_multiplier_fault;

        // Type-safe accessor — compile error if field not present
        template <SolutionFlags F>
        PYLITH_KERNEL const auto& get() const noexcept requires(Layout::has(F)) {
            if constexpr (F == SolutionFlags::INERTIA) { return velocity.member;}
            if constexpr (F == SolutionFlags::FAULT) { return lagrange_multiplier_fault.member;}
        } // get()

    }; // Unpacked

    /// Unpack solution fields from array into struct with names.
    template <typename Dim>
    PYLITH_KERNEL static Unpacked<Dim> unpack(const pylith::integer sOff[],
                                              const pylith::integer sOff_x[],
                                              const pylith::scalar s[],
                                              const pylith::scalar s_t[],
                                              const pylith::scalar s_x[]) {
        Unpacked<Dim> data;

        data.displacement = {
            &s[sOff[common::toIndex(Fields::DISPLACEMENT)]],
            s_t ? &s_t[sOff[common::toIndex(Fields::DISPLACEMENT)]] : nullptr,
            s_x ? &s_x[sOff_x[common::toIndex(Fields::DISPLACEMENT)]] : nullptr,
        };

        if constexpr (has(SolutionFlags::INERTIA)) {
            data.velocity = {
                &s[sOff[common::toIndex(Fields::INERTIA)]],
                s_t ? &s_t[sOff[common::toIndex(Fields::INERTIA)]] : nullptr,
                s_x ? &s_x[sOff_x[common::toIndex(Fields::INERTIA)]] : nullptr,
            };
        }

        if constexpr (has(SolutionFlags::FAULT)) {
            data.lagrange_multiplier_fault = {
                &s[sOff[common::toIndex(Fields::LAGRANGE_MULTIPLIER_FAULT)]],
                s_t ? &s_t[sOff[common::toIndex(Fields::LAGRANGE_MULTIPLIER_FAULT)]] : nullptr,
                s_x ? &s_x[sOff_x[common::toIndex(Fields::LAGRANGE_MULTIPLIER_FAULT)]] : nullptr,
            };
        }

        return data;
    } // unpack

}; // SolutionLayout
