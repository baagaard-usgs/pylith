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
#include "pylith/fekernels/common/Matrix.hh"
#include "pylith/fekernels/common/PetscJacobian.hh"

#include <cstddef>


namespace pylith::fekernels::momentum::stress::elasticity {
    template<size_t dim, class AuxiliaryLayout, class AuxiliaryUnpacked> struct IsotropicLinear;
} // namespace


/// Constitutive behavior for the isotropic linear elastic bulk rheology.
template<size_t dim, class AuxiliaryLayout, class AuxiliaryUnpacked>
struct pylith::fekernels::momentum::stress::elasticity::IsotropicLinear {
    /// σ = σ^mean + σ^dev.
    PYLITH_KERNEL static void cauchyStress(pylith::fekernels::common::Matrix<dim>& stress,
                                           const pylith::fekernels::common::Matrix<dim>& strain,
                                           const AuxiliaryUnpacked& auxiliary) {
        stress.zero();
        meanStress(stress, auxiliary.bulk_modulus(), strain);
        deviatoricStress(stress, auxiliary.shear_modulus(), strain);

        // Reference stress subtracted from residual
        if constexpr (AuxiliaryLayout::has(pylith::fekernels::pde::elasticity::ISOTROPIC_LINEAR_REFERENCE_STRESS)) {
            const auto& reference_stress = auxiliary.template get<pylith::fekernels::pde::elasticity::ISOTROPIC_LINEAR_REFERENCE_STRESS>();
            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {
                    stress(i, j) -= reference_stress(i,j);
                } // for
            } // for
        } // if
    } // compute

    /// J(f,g,df,dg) = C(f,df,g,dg)
    /// C_ijkl = λ δ_ij δ_kl + μ(δ_ik δ_jl + δ_il δ_jk)
    PYLITH_KERNEL static void cauchyStressTangent(pylith::scalar Jf3[],
                                                  const AuxiliaryUnpacked& auxiliary,
                                                  const pylith::scalar sign) noexcept {
        const pylith::scalar lambda = auxiliary.bulk_modulus() - 2.0/3.0 * auxiliary.shear_modulus();
        const pylith::scalar mu = auxiliary.shear_modulus();

        // C(f,df,g,dg)
        using Jacobian = pylith::fekernels::common::PetscJacobian<dim>;
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                for (size_t k = 0; k < dim; k++) {
                    for (size_t l = 0; l < dim; l++) {
                        Jf3[Jacobian::index3(i, j, k, l)] = sign * (
                            lambda * (i == j) * (k == l)
                            + mu * ((i == k)*(j == l) + (i == l)*(j == k)));
                    } // for
                } // for
            } // for
        } // for
    } // tangent

    /// σ^mean_ij = K * ε_kk = K * tr(ε)
    PYLITH_KERNEL static void meanStress(pylith::fekernels::common::Matrix<dim>& stress,
                                         const pylith::scalar bulkModulus,
                                         const pylith::fekernels::common::Matrix<dim>& strain) {
        assert(bulkModulus > 0.0);

        for (size_t i = 0; i < dim; ++i) {
            stress(i, i) += bulkModulus * strain.trace();
        } // for
    } // meanStress

    /// σ^dev_ij = 2 μ ε_ij - 1/3 tr(ε) δ_ij
    PYLITH_KERNEL static void deviatoricStress(pylith::fekernels::common::Matrix<dim>& stress,
                                               const pylith::scalar shearModulus,
                                               const pylith::fekernels::common::Matrix<dim>& strain) {
        assert(shearModulus > 0.0);

        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                stress(i, j) += 2.0*shearModulus * (strain(i,j) - (i == j)*1.0/3.0*strain.trace());
            } // for
        } // for
    } // meanStress

}; // IsotropicLinear
