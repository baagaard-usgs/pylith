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

#include <cstddef>
#include <map>
#include <stdexcept>
#include <tuple>
#include <string>

namespace pylith::fekernels::pde::elasticity::isotropic_linear {
    class DivStressRegistry;
} // namespace


class pylith::fekernels::pde::elasticity::isotropic_linear::DivStressRegistry {
private:

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

public:

    DivStressRegistry(void);
    inline
    ResidualFn f1(const size_t dim,
                  const StrainFlags strainFlags,
                  const SolutionFlags solutionFlags,
                  const MomentumFlags momentumFlags,
                  const AuxiliaryFlags auxiliaryFlags) const;

    inline
    JacobianFn Jf3uu(const size_t dim,
                     const StrainFlags strainFlags,
                     const SolutionFlags solutionFlags,
                     const MomentumFlags momentumFlags,
                     const AuxiliaryFlags auxiliaryFlags) const;

private:

    template<typename Dim, StrainFlags strainFlags, SolutionFlags solutionFlags, MomentumFlags momentumFlags, AuxiliaryFlags auxiliaryFlags>
    inline void _registerKernels(void);

private:

    ResidualFnMap _f1;
    JacobianFnMap _Jf3uu;

}; // DivStressRegistry


#include "DivStressRegistry.icc"
