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
#include "pylith/fekernels/common/Tensor2.hh"

#include <cstddef>


namespace pylith::fekernels::momentum::stress::elasticity {
    template<typename Dim, class AuxiliaryUnpacked> struct IsotropicLinear;
} // namespace


/// Constitutive behavior for the isotropic linear elastic bulk rheology.
template<typename Dim, class AuxiliaryLayout>
class pylith::fekernels::momentum::stress::elasticity::IsotropicLinear {
public:

    using AuxiliaryFlags = pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryFlags;
    using AuxiliaryUnpacked = typename AuxiliaryLayout::template Unpacked<Dim>;

    /// σ = σ^mean + σ^dev.
    PYLITH_KERNEL static void cauchyStress(pylith::fekernels::common::Tensor2<Dim>& stress,
                                           const pylith::fekernels::common::Tensor2<Dim>& strain,
                                           const AuxiliaryUnpacked& auxiliary) noexcept;

    /// J(f,g,df,dg) = C(f,df,g,dg)
    /// C_ijkl = λ δ_ij δ_kl + μ(δ_ik δ_jl + δ_il δ_jk)
    PYLITH_KERNEL static void cauchyStressTangent(pylith::scalar Jf3[],
                                                  const AuxiliaryUnpacked& auxiliary,
                                                  const pylith::scalar sign) noexcept;

    /// σ^mean_ij = K * ε_kk = K * tr(ε)
    PYLITH_KERNEL static void meanStress(pylith::fekernels::common::Tensor2<Dim>& stress,
                                         const pylith::scalar bulkModulus,
                                         const pylith::fekernels::common::Tensor2<Dim>& strain) noexcept;

    /// σ^dev_ij = 2 μ ε_ij - 1/3 tr(ε) δ_ij
    PYLITH_KERNEL static void deviatoricStress(pylith::fekernels::common::Tensor2<Dim>& stress,
                                               const pylith::scalar shearModulus,
                                               const pylith::fekernels::common::Tensor2<Dim>& strain) noexcept;

}; // IsotropicLinear
