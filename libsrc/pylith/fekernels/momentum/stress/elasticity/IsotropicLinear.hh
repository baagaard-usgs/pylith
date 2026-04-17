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
#include "pylith/fekernels/common/StorageTypes.hh"
#include "pylith/fekernels/common/PetscJacobian.hh"

#include <cstddef>


namespace pylith::fekernels::momentum::stress::elasticity {
    template<typename Dim, class AuxiliaryUnpacked> struct IsotropicLinear;
} // namespace


/// Constitutive behavior for the isotropic linear elastic bulk rheology.
template<typename Dim, class AuxiliaryLayout>
struct pylith::fekernels::momentum::stress::elasticity::IsotropicLinear {
    using AuxiliaryFlags = pylith::fekernels::pde::elasticity::isotropic_linear::AuxiliaryFlags;
    using AuxiliaryUnpacked = typename AuxiliaryLayout::template Unpacked<Dim>;

    /// σ = σ^mean + σ^dev.
    PYLITH_KERNEL static void cauchyStress(pylith::fekernels::common::Tensor2<Dim>& stress,
                                           const pylith::fekernels::common::Tensor2<Dim>& strain,
                                           const AuxiliaryUnpacked& auxiliary) {
        stress.zero();
        meanStress(stress, auxiliary.bulk_modulus(), strain);
        deviatoricStress(stress, auxiliary.shear_modulus(), strain);

        // Reference stress subtracted from residual
        if constexpr (requires { auxiliary.template get<AuxiliaryFlags::REFERENCE_STRESS>(); }) {
            const auto& reference_stress = auxiliary.template get<AuxiliaryFlags::REFERENCE_STRESS>();
            for (size_t i = 0; i < Dim::value; i++) {
                for (size_t j = 0; j < Dim::value; j++) {
                    stress(i, j) -= reference_stress(i,j);
                } // for
            } // for
        } // if
    } // cauchyStress

    /// J(f,g,df,dg) = C(f,df,g,dg)
    /// C_ijkl = λ δ_ij δ_kl + μ(δ_ik δ_jl + δ_il δ_jk)
    PYLITH_KERNEL static void cauchyStressTangent(pylith::scalar Jf3[],
                                                  const AuxiliaryUnpacked& auxiliary,
                                                  const pylith::scalar sign) noexcept {
        const pylith::scalar lambda = auxiliary.bulk_modulus() - 2.0/3.0 * auxiliary.shear_modulus();
        const pylith::scalar mu = auxiliary.shear_modulus();

        // C(f,df,g,dg)
        using Jacobian = pylith::fekernels::common::PetscJacobian<Dim>;
        for (size_t i = 0; i < Dim::value; i++) {
            for (size_t j = 0; j < Dim::value; j++) {
                for (size_t k = 0; k < Dim::value; k++) {
                    for (size_t l = 0; l < Dim::value; l++) {
                        Jf3[Jacobian::index3(i, j, k, l)] = sign * (
                            lambda * (i == j) * (k == l)
                            + mu * ((i == k)*(j == l) + (i == l)*(j == k)));
                    } // for
                } // for
            } // for
        } // for
    } // cauchyStressTangent

    /// σ^mean_ij = K * ε_kk = K * tr(ε)
    PYLITH_KERNEL static void meanStress(pylith::fekernels::common::Tensor2<Dim>& stress,
                                         const pylith::scalar bulkModulus,
                                         const pylith::fekernels::common::Tensor2<Dim>& strain) {
        assert(bulkModulus > 0.0);

        for (size_t i = 0; i < Dim::value; ++i) {
            stress(i, i) += bulkModulus * strain.trace();
        } // for
    } // meanStress

    /// σ^dev_ij = 2 μ ε_ij - 1/3 tr(ε) δ_ij
    PYLITH_KERNEL static void deviatoricStress(pylith::fekernels::common::Tensor2<Dim>& stress,
                                               const pylith::scalar shearModulus,
                                               const pylith::fekernels::common::Tensor2<Dim>& strain) {
        assert(shearModulus > 0.0);

        for (size_t i = 0; i < Dim::value; ++i) {
            for (size_t j = 0; j < Dim::value; ++j) {
                stress(i, j) += 2.0*shearModulus * (strain(i,j) - (i == j)*1.0/3.0*strain.trace());
            } // for
        } // for
    } // meanStress

}; // IsotropicLinear
