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
#include "pylith/fekernels/common/Flags.hh"
#include "pylith/fekernels/momentum/pde/MomentumLayout.hh"


#include <cassert>
#include <cstddef>


namespace pylith::fekernels::pde::elasticity {
    // Flags for optional subfields in isotropic linear auxiliary field
    enum IsotropicLinearFlags : size_t {
        ISOTROPIC_LINEAR_DEFAULT=0,
        ISOTROPIC_LINEAR_REFERENCE_STRESS=1 << 2,
        ISOTROPIC_LINEAR_REFERENCE_STRAIN=1 << 3,

    }; // IsotropicLinearFlags

    inline constexpr IsotropicLinearFlags
    operator|(IsotropicLinearFlags a,
              IsotropicLinearFlags b) noexcept {
        return static_cast<IsotropicLinearFlags>(
            static_cast<size_t>(a) | static_cast<size_t>(b));
    } // operator|


    inline constexpr IsotropicLinearFlags&
    operator|=(IsotropicLinearFlags& a,
               IsotropicLinearFlags b) noexcept {
        a = a | b;
        return a;
    }


    template<pylith::fekernels::momentum::MomentumFlags mflags, IsotropicLinearFlags flags> struct IsotropicLinearLayout;
} // namespace

/// Layout of auxiliary field for isotropic linear elasticity.
template<pylith::fekernels::momentum::MomentumFlags mflags,
         pylith::fekernels::pde::elasticity::IsotropicLinearFlags flags>
struct pylith::fekernels::pde::elasticity::IsotropicLinearLayout {
    // Is a given flag present?
    static constexpr bool has(pylith::fekernels::momentum::MomentumFlags f) {
        return (mflags & f);
    } // has

    // Is a given flag present?
    static constexpr bool has(IsotropicLinearFlags f) {
        return (flags & f);
    } // has

    // Optional subfields
    // Empty base class optimization ensures zero size when unused
    template <bool present, typename T>
    struct OptionalMember {};

    template <typename T>
    struct OptionalMember<true, T> { T member; };

    /// Order of auxiliary subfields.
#if 0 // Desired order but inconsistent with creation of auxiliary field
    enum Fields : int {
        DENSITY=0,
        BULK_MODULUS=1,
        SHEAR_MODULUS=2,
        BODY_FORCE=3,
        GRAVITY=3 + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0),
        REFERENCE_STRESS=3
                          + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                          + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0),
        REFERENCE_STRAIN=3
                          + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                          + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0)
                          + (has(ISOTROPIC_LINEAR_REFERENCE_STRESS) ? 1 : 0),
        NUM_FIELDS=3
                    + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                    + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0)
                    + (has(ISOTROPIC_LINEAR_REFERENCE_STRESS) ? 1 : 0)
                    +(has(ISOTROPIC_LINEAR_REFERENCE_STRAIN) ? 1 : 0),
#else
    enum Fields : int {
        DENSITY=0,
        BODY_FORCE=1,
        GRAVITY=1 + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0),
        REFERENCE_STRESS=1
                          + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                          + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0),
        REFERENCE_STRAIN=3
                          + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                          + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0)
                          + (has(ISOTROPIC_LINEAR_REFERENCE_STRESS) ? 1 : 0),
        BULK_MODULUS=1
                      + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                      + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0)
                      + (has(ISOTROPIC_LINEAR_REFERENCE_STRESS) ? 1 : 0)
                      +(has(ISOTROPIC_LINEAR_REFERENCE_STRAIN) ? 1 : 0),
        SHEAR_MODULUS=2
                       + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                       + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0)
                       +(has(ISOTROPIC_LINEAR_REFERENCE_STRESS) ? 1 : 0)
                       +(has(ISOTROPIC_LINEAR_REFERENCE_STRAIN) ? 1 : 0),
        NUM_FIELDS=2
                    + (has(pylith::fekernels::momentum::BODY_FORCE) ? 1 : 0)
                    + (has(pylith::fekernels::momentum::GRAVITY) ? 1 : 0)
                    + (has(ISOTROPIC_LINEAR_REFERENCE_STRESS) ? 1 : 0)
                    +(has(ISOTROPIC_LINEAR_REFERENCE_STRAIN) ? 1 : 0),
#endif
    }; // IsotropicLinearFlags

    /// Struct wth names for holding subfields
    template <size_t dim>
    struct Unpacked {
        using Layout = pylith::fekernels::pde::elasticity::IsotropicLinearLayout<mflags, flags>;

        pylith::fekernels::common::ScalarField density;
        pylith::fekernels::common::ScalarField bulk_modulus;
        pylith::fekernels::common::ScalarField shear_modulus;

        // Optional — zero size if absent
        [[no_unique_address]]
        OptionalMember<has(pylith::fekernels::momentum::BODY_FORCE), pylith::fekernels::common::VectorField<dim> > body_force;

        [[no_unique_address]]
        OptionalMember<has(pylith::fekernels::momentum::GRAVITY), pylith::fekernels::common::VectorField<dim> > gravitational_acceleration;

        [[no_unique_address]]
        OptionalMember<has(ISOTROPIC_LINEAR_REFERENCE_STRESS), pylith::fekernels::common::MatrixField<dim> > reference_stress;

        [[no_unique_address]]
        OptionalMember<has(ISOTROPIC_LINEAR_REFERENCE_STRAIN), pylith::fekernels::common::MatrixField<dim> > reference_strain;

        // Type-safe accessor — compile error if field not present
        template <pylith::fekernels::momentum::MomentumFlags F>
        auto& get() requires(Layout::has(F)) {
            if constexpr (F == pylith::fekernels::momentum::BODY_FORCE) { return body_force.member;}
            if constexpr (F == pylith::fekernels::momentum::GRAVITY) { return gravitational_acceleration.member;}
        } // get()

        template <IsotropicLinearFlags F>
        auto& get() requires(Layout::has(F)) {
            if constexpr (F == ISOTROPIC_LINEAR_REFERENCE_STRESS) { return reference_stress.member;}
            if constexpr (F == ISOTROPIC_LINEAR_REFERENCE_STRAIN) { return reference_strain.member;}
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

        data.density = {
            &s[sOff[DENSITY]],
            s_t ? &s_t[sOff[DENSITY]] : nullptr,
            s_x ? &s_x[sOff_x[DENSITY]] : nullptr,
        };

        data.bulk_modulus = {
            &s[sOff[BULK_MODULUS]],
            s_t ? &s_t[sOff[BULK_MODULUS]] : nullptr,
            s_x ? &s_x[sOff_x[BULK_MODULUS]] : nullptr,
        };

        data.shear_modulus = {
            &s[sOff[SHEAR_MODULUS]],
            s_t ? &s_t[sOff[SHEAR_MODULUS]] : nullptr,
            s_x ? &s_x[sOff_x[SHEAR_MODULUS]] : nullptr,
        };

        if constexpr (has(pylith::fekernels::momentum::BODY_FORCE)) {
            data.body_force = {
                &s[sOff[BODY_FORCE]],
                s_t ? &s_t[sOff[BODY_FORCE]] : nullptr,
                s_x ? &s_x[sOff_x[BODY_FORCE]] : nullptr,
            };
        }

        if constexpr (has(pylith::fekernels::momentum::GRAVITY)) {
            data.gravitational_acceleration = {
                &s[sOff[GRAVITY]],
                s_t ? &s_t[sOff[GRAVITY]] : nullptr,
                s_x ? &s_x[sOff_x[GRAVITY]] : nullptr,
            };
        }

        if constexpr (has(ISOTROPIC_LINEAR_REFERENCE_STRESS)) {
            data.reference_stress = {
                &s[sOff[REFERENCE_STRESS]],
                s_t ? &s_t[sOff[REFERENCE_STRESS]] : nullptr,
                s_x ? &s_x[sOff_x[REFERENCE_STRESS]] : nullptr,
            };
        }

        if constexpr (has(ISOTROPIC_LINEAR_REFERENCE_STRAIN)) {
            data.reference_strain = {
                &s[sOff[REFERENCE_STRAIN]],
                s_t ? &s_t[sOff[REFERENCE_STRAIN]] : nullptr,
                s_x ? &s_x[sOff_x[REFERENCE_STRAIN]] : nullptr,
            };
        }

        return data;
    } // unpack

}; // IsotropicLinearLayout
