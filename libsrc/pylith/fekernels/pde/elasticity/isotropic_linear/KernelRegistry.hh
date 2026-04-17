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

#include "pylith/fekernels/pde/elasticity/SolutionLayout.hh"
#include "pylith/fekernels/momentum/pde/MomentumLayout.hh"
#include "pylith/fekernels/momentum/pde/BodyForce.hh"
#include "pylith/fekernels/momentum/pde/DivStress.hh"
#include "pylith/fekernels/momentum/strain/StrainFlags.hh"
#include "pylith/fekernels/momentum/strain/InfinitesimalStrain.hh"
#include "pylith/fekernels/momentum/stress/elasticity/IsotropicLinear.hh"
#include "AuxiliaryLayout.hh"

#include <cassert>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <tuple>
#include <string>


namespace pylith::fekernels::pde::elasticity::isotropic_linear {
    class BodyForceRegistry;
    class DivStressRegistry;

} // namespace

class pylith::fekernels::pde::elasticity::isotropic_linear::BodyForceRegistry {
    typedef pylith::fekernels::momentum::MomentumFlags MomentumFlags;
    typedef pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryFlags AuxiliaryFlags;

    template<MomentumFlags mflags, AuxiliaryFlags aflags>
    using AuxiliaryLayout = pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryLayout<mflags, aflags>;

    using ResidualFn = PetscPointFn*;
    using ResidualFnMap = std::map<std::tuple<size_t, MomentumFlags, AuxiliaryFlags>, ResidualFn>;

    ResidualFnMap _f0;

    template<typename Dim, MomentumFlags momentumFlags, AuxiliaryFlags auxiliaryFlags>
    inline void _registerKernels(void) {
        _f0[std::make_tuple(Dim::value, momentumFlags, auxiliaryFlags)] =
            pylith::fekernels::momentum::BodyForce< Dim, AuxiliaryLayout<momentumFlags, auxiliaryFlags> >::f0;
    }

public:

    BodyForceRegistry(void) {
        // 2D
        _registerKernels<Dim2, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();

        _registerKernels<Dim2, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS>();
        _registerKernels<Dim2, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS>();
        _registerKernels<Dim2, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS>();
        _registerKernels<Dim2, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS>();

        // 3D
        _registerKernels<Dim3, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();

        _registerKernels<Dim3, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS>();
        _registerKernels<Dim3, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS>();
        _registerKernels<Dim3, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS>();
        _registerKernels<Dim3, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS>();

    } // constructor

    inline
    ResidualFn f0(const size_t dim,
                  const MomentumFlags momentumFlags,
                  const AuxiliaryFlags auxiliaryFlags) const {
        auto iter = _f0.find(std::tuple(dim, momentumFlags, auxiliaryFlags));
        if (iter == _f0.end()) {
            throw std::invalid_argument(unsupportedMessage(dim, momentumFlags, auxiliaryFlags));
        } // if
        return iter->second;
    } // get

    static inline
    std::string unsupportedMessage(const size_t dim,
                                   const MomentumFlags momentumFlags,
                                   const AuxiliaryFlags auxiliaryFlags) {
        return "Unsupported physics combination: "
               " dim=" + std::to_string(dim)
               + ", momentum flags=" + std::to_string(uint32_t(momentumFlags))
               + ", auxiliary flags=" + std::to_string(uint32_t(auxiliaryFlags));
    } // unsupportedMessage

}; // BodyForceRegistry


class pylith::fekernels::pde::elasticity::isotropic_linear::DivStressRegistry {
    typedef pylith::fekernels::pde::elasticity::SolutionFlags SolutionFlags;
    typedef pylith::fekernels::momentum::StrainFlags StrainFlags;
    typedef pylith::fekernels::momentum::MomentumFlags MomentumFlags;
    typedef pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryFlags AuxiliaryFlags;

    template<MomentumFlags mflags, AuxiliaryFlags aflags>
    using AuxiliaryLayout = pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryLayout<mflags, aflags>;

    template<SolutionFlags flags>
    using SolutionLayout = pylith::fekernels::pde::elasticity::SolutionLayout<flags>;

    // template<typename Dim, class StrainModel, class StressModel, class SolutionLayout, class AuxiliaryLayout>


    using ResidualFn = PetscPointFn*;
    using ResidualFnMap = std::map<std::tuple<size_t, StrainFlags, SolutionFlags, MomentumFlags, AuxiliaryFlags>, ResidualFn>;
    using JacobianFn = PetscPointJacFn*;
    using JacobianFnMap = std::map<std::tuple<size_t, StrainFlags, SolutionFlags, MomentumFlags, AuxiliaryFlags>, JacobianFn>;

    ResidualFnMap _f1;
    JacobianFnMap _Jf3uu;

    template<typename Dim, StrainFlags strainFlags, SolutionFlags solutionFlags, MomentumFlags momentumFlags, AuxiliaryFlags auxiliaryFlags>
    inline void _registerKernels(void) {
        using SolutionLayout = pylith::fekernels::pde::elasticity::SolutionLayout<solutionFlags>;
        using StrainModel = pylith::fekernels::momentum::InfinitesimalStrain<Dim, SolutionLayout>;
        using StressModel = pylith::fekernels::momentum::stress::elasticity::IsotropicLinear< Dim, AuxiliaryLayout<momentumFlags, auxiliaryFlags> >;
        _f1[std::make_tuple(Dim::value, strainFlags, solutionFlags, momentumFlags, auxiliaryFlags)] =
            pylith::fekernels::momentum::DivStress<Dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout<momentumFlags, auxiliaryFlags> >::f1;

        _Jf3uu[std::make_tuple(Dim::value, strainFlags, solutionFlags, momentumFlags, auxiliaryFlags)] =
            pylith::fekernels::momentum::DivStress<Dim, StrainModel, StressModel, SolutionLayout, AuxiliaryLayout<momentumFlags, auxiliaryFlags> >::Jf3uu;
    } // registerKernels

public:

    DivStressRegistry(void) {
        // 2D
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();

        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();

        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim2, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();

        // 3D
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::DEFAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();

        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();

        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::DEFAULT, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::DEFAULT, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::DEFAULT>();
        _registerKernels<Dim3, StrainFlags::INFINITESIMAL, SolutionFlags::FAULT | SolutionFlags::INERTIA, MomentumFlags::BODY_FORCE | MomentumFlags::GRAVITY, AuxiliaryFlags::REFERENCE_STRESS|AuxiliaryFlags::REFERENCE_STRAIN>();

    } // constructor

    inline
    ResidualFn f1(const size_t dim,
                  const StrainFlags strainFlags,
                  const SolutionFlags solutionFlags,
                  const MomentumFlags momentumFlags,
                  const AuxiliaryFlags auxiliaryFlags) const {
        auto iter = _f1.find(std::tuple(dim, strainFlags, solutionFlags, momentumFlags, auxiliaryFlags));
        if (iter == _f1.end()) {
            throw std::invalid_argument(unsupportedMessage(dim, strainFlags, solutionFlags, momentumFlags, auxiliaryFlags));
        } // if
        return iter->second;
    } // get

    inline
    JacobianFn Jf3uu(const size_t dim,
                     const StrainFlags strainFlags,
                     const SolutionFlags solutionFlags,
                     const MomentumFlags momentumFlags,
                     const AuxiliaryFlags auxiliaryFlags) const {
        auto iter = _Jf3uu.find(std::tuple(dim, strainFlags, solutionFlags, momentumFlags, auxiliaryFlags));
        if (iter == _Jf3uu.end()) {
            throw std::invalid_argument(unsupportedMessage(dim, strainFlags, solutionFlags, momentumFlags, auxiliaryFlags));
        } // if
        return iter->second;
    } // get

    static inline
    std::string unsupportedMessage(const size_t dim,
                                   const StrainFlags strainFlags,
                                   const SolutionFlags solutionFlags,
                                   const MomentumFlags momentumFlags,
                                   const AuxiliaryFlags auxiliaryFlags) {
        return "Unsupported physics combination: "
               " dim=" + std::to_string(dim)
               + ", strain flags=" + std::to_string(uint32_t(strainFlags))
               + ", solution flags=" + std::to_string(uint32_t(solutionFlags))
               + ", momentum flags=" + std::to_string(uint32_t(momentumFlags))
               + ", auxiliary flags=" + std::to_string(uint32_t(auxiliaryFlags));
    } // unsupportedMessage

}; // DivStressRegistry
