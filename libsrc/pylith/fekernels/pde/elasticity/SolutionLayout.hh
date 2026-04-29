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

#include "pylith/fekernels/common/Kernel.hh"
#include "pylith/fekernels/common/Utils.hh"
#include "pylith/fekernels/common/VectorField.hh"
#include "pylith/fekernels/common/OptionalFields.hh"

#include <cstddef>


namespace pylith::fekernels::pde::elasticity {
    // Flags for elasticity solution
    enum class SolutionFlags : size_t {
        DEFAULT=0,
        FAULT=1 << 0,
        INERTIA=1 << 1,
    }; // SolutionFlags

    PYLITH_KERNEL constexpr SolutionFlags operator|(SolutionFlags a,
                                                    SolutionFlags b) noexcept;


    PYLITH_KERNEL constexpr SolutionFlags&operator|=(SolutionFlags& a,
                                                     SolutionFlags b) noexcept;


    PYLITH_KERNEL constexpr bool operator&(SolutionFlags a,
                                           SolutionFlags b) noexcept;

    template<SolutionFlags flags> class SolutionLayout;
} // namespace


/// Layout of solution field for elasticity
template<pylith::fekernels::pde::elasticity::SolutionFlags flags>
class pylith::fekernels::pde::elasticity::SolutionLayout {
public:

    /// Order of solution subfields
    enum class Fields : uint32_t {
        DISPLACEMENT=0,
        INERTIA=1,
        LAGRANGE_MULTIPLIER_FAULT=1 + pylith::fekernels::common::addIf(has(SolutionFlags::INERTIA)),
        NUM_FIELDS=1 + pylith::fekernels::common::addIf(has(SolutionFlags::INERTIA))
                    + pylith::fekernels::common::addIf(has(SolutionFlags::FAULT))
    };

    // Is a given flag present?
    static constexpr bool has(SolutionFlags f) noexcept;

    /// Struct wth names for holding subfields
    template <typename Dim>
    class Unpacked {
public:

        using Layout = pylith::fekernels::pde::elasticity::SolutionLayout<flags>;

        pylith::fekernels::common::VectorField<Dim> displacement;

        // Optional — zero size if absent
        [[no_unique_address]]
        pylith::fekernels::common::OptionalMember<has(SolutionFlags::INERTIA), pylith::fekernels::common::VectorField<Dim> > velocity;

        [[no_unique_address]]
        pylith::fekernels::common::OptionalMember<has(SolutionFlags::FAULT), pylith::fekernels::common::VectorField<Dim> > lagrange_multiplier_fault;

        // Type-safe accessor — compile error if field not present
        template <SolutionFlags F>
        PYLITH_KERNEL const auto& get() const noexcept requires(Layout::has(F));

    }; // Unpacked

    /// Unpack solution fields from array into struct with names.
    template <typename Dim>
    PYLITH_KERNEL static Unpacked<Dim> unpack(const pylith::integer sOff[],
                                              const pylith::integer sOff_x[],
                                              const pylith::scalar s[],
                                              const pylith::scalar s_t[],
                                              const pylith::scalar s_x[]);

}; // SolutionLayout


#include "SolutionLayout.icc"
