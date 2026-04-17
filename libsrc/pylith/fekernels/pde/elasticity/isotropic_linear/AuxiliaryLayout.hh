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
#include "pylith/fekernels/momentum/pde/MomentumLayout.hh"


#include <cassert>
#include <cstddef>

namespace common = pylith::fekernels::common;

namespace pylith::fekernels::pde::elasticity::isotropic_linear {
    // Flags for optional subfields in isotropic linear auxiliary field
    enum class AuxiliaryFlags : size_t {
        DEFAULT=0,
        REFERENCE_STRESS=1 << 2,
        REFERENCE_STRAIN=1 << 3,

    }; // Flags

    PYLITH_KERNEL constexpr AuxiliaryFlags
    operator|(AuxiliaryFlags a,
              AuxiliaryFlags b) noexcept {
        return static_cast<AuxiliaryFlags>(
            static_cast<size_t>(a) | static_cast<size_t>(b));
    } // operator|


    PYLITH_KERNEL constexpr AuxiliaryFlags&
    operator|=(AuxiliaryFlags& a,
               AuxiliaryFlags b) noexcept {
        a = a | b;
        return a;
    }


    PYLITH_KERNEL constexpr bool
    operator&(AuxiliaryFlags a,
              AuxiliaryFlags b) noexcept {
        return static_cast<size_t>(a) & static_cast<size_t>(b);
    }


    template<pylith::fekernels::momentum::MomentumFlags mflags, AuxiliaryFlags flags> struct AuxiliaryLayout;
} // namespace


/// Layout of auxiliary field for isotropic linear elasticity.
template<pylith::fekernels::momentum::MomentumFlags mflags,
         pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryFlags flags>
struct pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryLayout {
    using MomentumFlags = pylith::fekernels::momentum::MomentumFlags;

    // Is a given flag present?
    PYLITH_KERNEL static constexpr bool has(MomentumFlags f) {
        return (mflags & f);
    } // has

    // Is a given flag present?
    PYLITH_KERNEL static constexpr bool has(AuxiliaryFlags f) {
        return (flags & f);
    } // has

    /// Order of auxiliary subfields.
#if 0 // Desired order but inconsistent with creation of auxiliary field
    enum Fields : int {
        DENSITY=0,
        BULK_MODULUS=1,
        SHEAR_MODULUS=2,
        BODY_FORCE=3,
        GRAVITY=3 + (has(MomentumFlags::BODY_FORCE) ? 1 : 0),
        AuxiliaryFlags::REFERENCE_STRESS=3
                                          + (has(pylith::fekernels::momentum::MomentumFlags::BODY_FORCE) ? 1 : 0)
                                          + (has(pylith::fekernels::momentum::MomentumFlags::GRAVITY) ? 1 : 0),
        REFERENCE_STRAIN=3
                          + (has(pylith::fekernels::momentum::MomentumFlags::BODY_FORCE) ? 1 : 0)
                          + (has(pylith::fekernels::momentum::MomentumFlags::GRAVITY) ? 1 : 0)
                          + (has(AuxiliaryFlags::REFERENCE_STRESS) ? 1 : 0),
        NUM_FIELDS=3
                    + (has(pylith::fekernels::momentum::MomentumFlags::BODY_FORCE) ? 1 : 0)
                    + (has(pylith::fekernels::momentum::MomentumFlags::GRAVITY) ? 1 : 0)
                    + (has(AuxiliaryFlags::REFERENCE_STRESS) ? 1 : 0)
                    +(has(AuxiliaryFlags::REFERENCE_STRAIN) ? 1 : 0),
    }; // Fields

#else
    enum class Fields : uint32_t {
        DENSITY=0,
        BODY_FORCE=1,
        GRAVITY=1 + common::addIf(has(MomentumFlags::BODY_FORCE)),
        REFERENCE_STRESS=1
                          + common::addIf(has(MomentumFlags::BODY_FORCE))
                          + common::addIf(has(MomentumFlags::GRAVITY)),
        REFERENCE_STRAIN=1
                          + common::addIf(has(MomentumFlags::BODY_FORCE))
                          + common::addIf(has(MomentumFlags::GRAVITY))
                          + common::addIf(has(AuxiliaryFlags::REFERENCE_STRESS)),
        BULK_MODULUS=1
                      + common::addIf(has(MomentumFlags::BODY_FORCE))
                      + common::addIf(has(MomentumFlags::GRAVITY))
                      + common::addIf(has(AuxiliaryFlags::REFERENCE_STRESS))
                      + common::addIf(has(AuxiliaryFlags::REFERENCE_STRAIN)),
        SHEAR_MODULUS=2
                       + common::addIf(has(MomentumFlags::BODY_FORCE))
                       + common::addIf(has(MomentumFlags::GRAVITY))
                       + common::addIf(has(AuxiliaryFlags::REFERENCE_STRESS))
                       + common::addIf(has(AuxiliaryFlags::REFERENCE_STRAIN)),
        NUM_FIELDS=3
                    + common::addIf(has(MomentumFlags::BODY_FORCE))
                    + common::addIf(has(MomentumFlags::GRAVITY))
                    + common::addIf(has(AuxiliaryFlags::REFERENCE_STRESS))
                    + common::addIf(has(AuxiliaryFlags::REFERENCE_STRAIN)),
    }; // Fields

#endif

    /// Struct wth names for holding subfields
    template <typename Dim>
    struct Unpacked {
        using Layout = pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryLayout<mflags, flags>;

        common::ScalarField density;
        common::ScalarField bulk_modulus;
        common::ScalarField shear_modulus;

        // Optional — zero size if absent
        [[no_unique_address]]
        common::OptionalMember<has(MomentumFlags::BODY_FORCE), common::VectorField<Dim> > body_force;

        [[no_unique_address]]
        common::OptionalMember<has(MomentumFlags::GRAVITY), common::VectorField<Dim> > gravitational_acceleration;

        [[no_unique_address]]
        common::OptionalMember<has(AuxiliaryFlags::REFERENCE_STRESS), common::Tensor2Field<Dim> > reference_stress;

        [[no_unique_address]]
        common::OptionalMember<has(AuxiliaryFlags::REFERENCE_STRAIN), common::Tensor2Field<Dim> > reference_strain;

        // Type-safe accessor — compile error if field not present
        template <MomentumFlags F>
        PYLITH_KERNEL auto& get() const noexcept requires(Layout::has(F)) {
            if constexpr (F == MomentumFlags::BODY_FORCE) { return body_force.member;}
            if constexpr (F == MomentumFlags::GRAVITY) { return gravitational_acceleration.member;}
        } // get()

        template <AuxiliaryFlags F>
        PYLITH_KERNEL const auto& get() const noexcept requires(Layout::has(F)) {
            if constexpr (F == AuxiliaryFlags::REFERENCE_STRESS) { return reference_stress.member;}
            if constexpr (F == AuxiliaryFlags::REFERENCE_STRAIN) { return reference_strain.member;}
        } // get()

    }; // Unpacked

    /// Unpack solution fields from array into struct with names.
    template <typename Dim>
    PYLITH_KERNEL static Unpacked<Dim>unpack(const pylith::integer sOff[],
                                             const pylith::integer sOff_x[],
                                             const pylith::scalar s[],
                                             const pylith::scalar s_t[],
                                             const pylith::scalar s_x[]) {
        Unpacked<Dim> data;

        data.density = {
            &s[sOff[common::toIndex(Fields::DENSITY)]],
            s_t ? &s_t[sOff[common::toIndex(Fields::DENSITY)]] : nullptr,
            s_x ? &s_x[sOff_x[common::toIndex(Fields::DENSITY)]] : nullptr,
        };

        data.bulk_modulus = {
            &s[sOff[common::toIndex(Fields::BULK_MODULUS)]],
            s_t ? &s_t[sOff[common::toIndex(Fields::BULK_MODULUS)]] : nullptr,
            s_x ? &s_x[sOff_x[common::toIndex(Fields::BULK_MODULUS)]] : nullptr,
        };

        data.shear_modulus = {
            &s[sOff[common::toIndex(Fields::SHEAR_MODULUS)]],
            s_t ? &s_t[sOff[common::toIndex(Fields::SHEAR_MODULUS)]] : nullptr,
            s_x ? &s_x[sOff_x[common::toIndex(Fields::SHEAR_MODULUS)]] : nullptr,
        };

        if constexpr (has(pylith::fekernels::momentum::MomentumFlags::BODY_FORCE)) {
            data.body_force = {
                &s[sOff[common::toIndex(Fields::BODY_FORCE)]],
                s_t ? &s_t[sOff[common::toIndex(Fields::BODY_FORCE)]] : nullptr,
                s_x ? &s_x[sOff_x[common::toIndex(Fields::BODY_FORCE)]] : nullptr,
            };
        }

        if constexpr (has(pylith::fekernels::momentum::MomentumFlags::GRAVITY)) {
            data.gravitational_acceleration = {
                &s[sOff[common::toIndex(Fields::GRAVITY)]],
                s_t ? &s_t[sOff[common::toIndex(Fields::GRAVITY)]] : nullptr,
                s_x ? &s_x[sOff_x[common::toIndex(Fields::GRAVITY)]] : nullptr,
            };
        }

        if constexpr (has(AuxiliaryFlags::REFERENCE_STRESS)) {
            data.reference_stress = {
                &s[sOff[common::toIndex(Fields::REFERENCE_STRESS)]],
                s_t ? &s_t[sOff[common::toIndex(Fields::REFERENCE_STRESS)]] : nullptr,
                s_x ? &s_x[sOff_x[common::toIndex(Fields::REFERENCE_STRESS)]] : nullptr,
            };
        }

        if constexpr (has(AuxiliaryFlags::REFERENCE_STRAIN)) {
            data.reference_strain = {
                &s[sOff[common::toIndex(Fields::REFERENCE_STRAIN)]],
                s_t ? &s_t[sOff[common::toIndex(Fields::REFERENCE_STRAIN)]] : nullptr,
                s_x ? &s_x[sOff_x[common::toIndex(Fields::REFERENCE_STRAIN)]] : nullptr,
            };
        }

        return data;
    } // unpack

}; // IsotropicLinearLayout
