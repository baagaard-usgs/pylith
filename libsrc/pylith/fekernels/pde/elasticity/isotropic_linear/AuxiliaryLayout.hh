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
#include "pylith/fekernels/common/ScalarField.hh"
#include "pylith/fekernels/common/VectorField.hh"
#include "pylith/fekernels/common/Tensor2Field.hh"
#include "pylith/fekernels/common/OptionalFields.hh"
#include "pylith/fekernels/momentum/pde/MomentumLayout.hh"


#include <cstddef>

namespace common = pylith::fekernels::common;

namespace pylith::fekernels::pde::elasticity::isotropic_linear {
    // Flags for optional subfields in isotropic linear auxiliary field
    enum class AuxiliaryFlags : size_t {
        DEFAULT=0,
        REFERENCE_STRESS=1 << 0,
        REFERENCE_STRAIN=1 << 1,

    }; // Flags

    PYLITH_KERNEL constexpr AuxiliaryFlags operator|(AuxiliaryFlags a,
                                                     AuxiliaryFlags b) noexcept;

    PYLITH_KERNEL constexpr AuxiliaryFlags&operator|=(AuxiliaryFlags& a,
                                                      AuxiliaryFlags b) noexcept;

    PYLITH_KERNEL constexpr bool operator&(AuxiliaryFlags a,
                                           AuxiliaryFlags b) noexcept;


    template<pylith::fekernels::momentum::MomentumFlags mflags, AuxiliaryFlags flags> class AuxiliaryLayout;
} // namespace


/// Layout of auxiliary field for isotropic linear elasticity.
template<pylith::fekernels::momentum::MomentumFlags mflags,
         pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryFlags flags>
class pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryLayout {
public:

    using MomentumFlags = pylith::fekernels::momentum::MomentumFlags;

    // Is a given flag present?
    PYLITH_KERNEL static constexpr bool has(MomentumFlags f);

    // Is a given flag present?
    PYLITH_KERNEL static constexpr bool has(AuxiliaryFlags f);

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
    class Unpacked {
public:

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
        PYLITH_KERNEL auto& get() const noexcept requires(Layout::has(F));

        template <AuxiliaryFlags F>
        PYLITH_KERNEL const auto& get() const noexcept requires(Layout::has(F));

    }; // Unpacked

    /// Unpack solution fields from array into struct with names.
    template <typename Dim>
    PYLITH_KERNEL static Unpacked<Dim>unpack(const pylith::integer sOff[],
                                             const pylith::integer sOff_x[],
                                             const pylith::scalar s[],
                                             const pylith::scalar s_t[],
                                             const pylith::scalar s_x[]);

}; // AuxiliaryLayout


#include "AuxiliaryLayout.icc"
