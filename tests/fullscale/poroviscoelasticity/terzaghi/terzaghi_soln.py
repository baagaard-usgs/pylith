# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
# @file tests/fullscale/poroelasticity/terzaghi/terzaghi_soln.py
#
# @brief Analytical solution to Terzaghi's problem.
#
# This is based on Cheng (Poroelasticity, Section 7.3.1) and its implementation in PETSc tutorial
# src/ts/tutorials/ex53.c. It is adapted for the poro-viscoelastic case using the correspondence 
# principle.
#
# Notes:
# - The traction loading is impulsive, so we use a custom PETSc TS (pylith/utils/TSAdaptImpulse),
#   which uses a very small time step for the first step before using the user-specified time step.
# - The accuracy of the solution is poor in the first few time steps due to the impulsive loading
#   so we only check the last time step in the simulation, which is more accurate and select a
#   a tolerance appropriate for the discretization size. 
#
#          Uy=0
#       ----------
#       |        |
# Ux=0  |        | Ux=0
#       |        |
#       |        |
#       ----------
#         +1 Pa
#
# Dirichlet boundary conditions
#   Ux(+-1,0) = 0
#   Uy(x,+1) = 0
# Neumann boundary conditions
#   \tau_normal(x,+1) = 1*Pa

import numpy
import mpmath

# Physical properties
rho_s = 2500  # kg / m**3
rho_f = 1000  # kg / m**3
mu_f = 1.0  # Pa*s
mu_s = 1.0e+18 # Pa*s
G = 1.0e+9  # Pa
K_sg = 10.0e+9  # Pa
K_fl = 2.0e+9  # Pa
K_d = 10.0e+9  # Pa
alpha = 0.8
phi = 0.1
k = 1.5e-13  # m**2

y_max = 10.0e+3  # m
y_min = 0.0  # m
P_0 = 1.0e+4  # Pa

# Height of column, m
H = y_max - y_min

M = 1.0 / (phi / K_fl + (alpha - phi) / K_sg)  # Pa
K_u = K_d + alpha * alpha * M  # Pa,      Cheng (B.5)
nu = (3.0 * K_d - 2.0 * G) / (2.0 * (3.0 * K_d + G))  # -,       Cheng (B.8)
nu_u = (3.0 * K_u - 2.0 * G) / (2.0 * (3.0 * K_u + G))  # -,       Cheng (B.9)
eta = (3.0 * alpha * G) / (3.0 * K_d + 4.0 * G)  # -,       Cheng (B.11)
S = (3.0 * K_u + 4.0 * G) / (M * (3.0 * K_d + 4.0 * G))  # Pa^{-1}, Cheng (B.14)
c = (k / mu_f) / S  # m^2 / s, Cheng (B.16)

dt0 = 1.0e-6 * 1.0e+8
dt = 5.0e+4
t_end = 100.0e+4
time = numpy.concatenate((numpy.array([dt0]), dt0+numpy.arange(dt, t_end+0.5*dt, dt)))
tsteps = numpy.array([time[-1]])

# ----------------------------------------------------------------------


class AnalyticalSoln(object):
    """Analytical solution to Terzaghi's problem"""

    SPACE_DIM = 2
    TENSOR_SIZE = 4
    ITERATIONS = 1600

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "pressure": self.pressure,
            "trace_strain": self.trace_strain,
            "porosity": self.porosity,
            "solid_density": self.solid_density,
            "fluid_density": self.fluid_density,
            "fluid_viscosity": self.fluid_viscosity,
            "solid_viscosity": self.solid_viscosity,
            "shear_modulus": self.shear_modulus,
            "drained_bulk_modulus": self.drained_bulk_modulus,
            "biot_coefficient": self.biot_coefficient,
            "biot_modulus": self.biot_modulus,
            "isotropic_permeability": self.isotropic_permeability,
            "initial_amplitude": {
                "bc_xneg": self.zero_vector,
                "bc_xpos": self.zero_vector,
                "bc_ypos_traction": self.bc_ypos_traction,
                "bc_ypos_pressure": self.zero_scalar,
                "bc_yneg": self.zero_vector,
            },
        }
        return

    def getField(self, name, mesh_entity, pts):
        if name in "initial_amplitude":
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def zero_scalar(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, 1), dtype=numpy.float64)

    def zero_vector(self, locs):
        (npts, dim) = locs.shape
        return numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)

    def solid_density(self, locs):
        """Compute solid_density field at locations."""
        (npts, dim) = locs.shape
        solid_density = rho_s * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return solid_density

    def fluid_density(self, locs):
        """Compute fluid density field at locations."""
        (npts, dim) = locs.shape
        fluid_density = rho_f * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_density

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations."""
        (npts, dim) = locs.shape
        shear_modulus = G * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def porosity(self, locs):
        """Compute porosity field at locations."""
        (npts, dim) = locs.shape
        porosity = phi * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return porosity

    def fluid_viscosity(self, locs):
        """Compute fluid_viscosity field at locations."""
        (npts, dim) = locs.shape
        fluid_viscosity = mu_f * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return fluid_viscosity

    def solid_viscosity(self, locs):
        """Compute solid_viscosity field at locations."""
        (npts, dim) = locs.shape
        solid_viscosity = mu_s * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return solid_viscosity

    def drained_bulk_modulus(self, locs):
        """Compute undrained bulk modulus field at locations."""
        (npts, dim) = locs.shape
        undrained_bulk_modulus = K_d * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return undrained_bulk_modulus

    def biot_coefficient(self, locs):
        """Compute biot coefficient field at locations."""
        (npts, dim) = locs.shape
        biot_coefficient = alpha * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_coefficient

    def biot_modulus(self, locs):
        """Compute biot modulus field at locations."""
        (npts, dim) = locs.shape
        biot_modulus = M * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return biot_modulus

    def isotropic_permeability(self, locs):
        """Compute isotropic permeability field at locations."""
        (npts, dim) = locs.shape
        isotropic_permeability = k * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return isotropic_permeability

    def displacement(self, locs):
        """Compute displacement field at locations."""
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        displacement = numpy.zeros((ntpts, npts, dim), dtype=numpy.float64)
        y_star = 1.0 - locs[:, 1] / H

        # Laplace transforms of displacement
        def u_t0_LT(y_star, s):
            G_bar = G * s / (s + (G / mu_s))
            K_bar = K_d * s / (s + (G / mu_s))
            alpha_bar = alpha * s / (s + (G / mu_s))

            inverse_M_bar = (alpha_bar - phi) / K_sg + phi / K_fl
            Ku_bar = K_bar + alpha_bar**2 * (1. / inverse_M_bar)
            nu_u_bar = (3.0 * Ku_bar - 2.0 * G_bar)/(2 * (3.0 * Ku_bar + G_bar))

            return (-(P_0 * H * (1.0 - 2.0 * nu_u_bar )) / (s * 2 * G_bar * (1.0 - nu_u_bar))) * (1.0 - y_star)

        def u_LT(y_star, s):
            G_bar = G * s / (s + (G / mu_s))
            K_bar = K_d * s / (s + (G / mu_s))
            alpha_bar = alpha * s / (s + (G / mu_s))

            inverse_M_bar = (alpha_bar - phi) / K_sg + phi / K_fl
            Ku_bar = K_bar + alpha_bar**2 * (1. / inverse_M_bar)
            nu_bar = (3.0 * K_bar - 2.0 * G_bar) / (2.0 * (3.0 * K_bar + G_bar))
            nu_u_bar = (3.0 * Ku_bar - 2.0 * G_bar)/(2 * (3.0 * Ku_bar + G_bar))

            return (-(P_0 * H * (1.0 - 2.0 * nu_u_bar )) / (s * 2 * G_bar * (1.0 - nu_u_bar))) * (1.0 - y_star) + (( P_0 *H *(nu_u_bar - nu_bar)) / (2.0 * G_bar * (1.0 - nu_u_bar) * (1.0 - nu_bar))) * self.F2_LT(y_star, s)

        for i_y, y in enumerate(y_star):
            for i_t, t in enumerate(tsteps):
                if t < 0.0:
                    displacement[0, i_y, 1] = mpmath.invertlaplace(lambda s: u_t0_LT(y, s), 0.0)
                else:
                    t_star = (c * t) / ((2 * H) ** 2)
                    displacement[i_t, i_y, 1] = mpmath.invertlaplace(lambda s: u_LT(y, s), t_star)

        return displacement

    def pressure(self, locs):
        """Compute pressure field at locations."""
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        pressure = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        y_star = 1.0 - locs[:, 1] / H
        mpmath.dps = 30; mpmath.pretty = True

        # Laplace transform of pressure
        def p_LT(y_star, s):
            G_bar = G * s / (s + (G / mu_s))
            K_bar = K_d * s / (s + (G / mu_s))
            alpha_bar = alpha * s / (s + (G / mu_s))

            eta_bar = (3.0 * alpha_bar * G_bar) / (3.0 * K_bar + 4.0 * G_bar)
            inverse_M_bar = (alpha_bar - phi) / K_sg + phi / K_fl
            S_bar = inverse_M_bar + 3.0 * alpha_bar**2 / (3.0 * K_bar + 4.0 * G_bar)

            return ((P_0 * eta_bar) /(G_bar * S_bar)) * self.F1_LT(y_star, s)

        for i_y, y in enumerate(y_star):
            for i_t, t in enumerate(tsteps):
                t_star = (c * t) / (4.0 * H**2)
                pressure[i_t, i_y, 0] = mpmath.invertlaplace(lambda s: p_LT(y, s), t_star)

        return pressure

    def trace_strain(self, locs):
        """Compute trace strain field at locations."""
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        trace_strain = numpy.zeros((ntpts, npts, 1), dtype=numpy.float64)
        y_star = 1.0 - locs[:, 1] / H

        # Laplace transform of trace strain
        def strain_LT(y_star, s):
            G_bar = G * s / (s + (G / mu_s))
            K_bar = K_d * s / (s + (G / mu_s))
            alpha_bar = alpha * s / (s + (G / mu_s))

            inverse_M_bar = (alpha_bar - phi) / K_sg + phi / K_fl
            Ku_bar = K_bar + alpha_bar**2 * (1. / inverse_M_bar)
            nu_bar = (3.0 * K_bar - 2.0 * G_bar) / (2.0 * (3.0 * K_bar + G_bar))
            nu_u_bar = (3.0 * Ku_bar - 2.0 * G_bar)/(2 * (3.0 * Ku_bar + G_bar))

            return -((P_0 * H * (1-2*nu_u_bar))/(2.0* s * G_bar * (1.0 - nu_u_bar) * H)) + (( P_0 * H * (nu_u_bar - nu_bar)) / (2.0 * G_bar * (1.0 - nu_u_bar) * (1.0 - nu_bar))) * self.F3_LT(y_star, s)

        for i_y, y in enumerate(y_star):
            for i_t, t in enumerate(tsteps):
                t_star = (c * t) / (4 * H**2)
                trace_strain[i_t, i_y, 0] = mpmath.invertlaplace(lambda s: strain_LT(y, s), t_star)

        return trace_strain

    # Series functions

    def F1_LT(self, y_star, s):
        # m = numpy.arange(1, 2 * self.ITERATIONS + 1, 2)
        # F1 = numpy.zeros(y_star.shape)
        # for i_y, y in enumerate(y_star):
        #     F1[i_y] = numpy.sum(4.0 / (m * numpy.pi) * mpmath.sin(0.5 * m * numpy.pi * y) * (1.0 / (s + (m * numpy.pi)**2)))
        # return F1

        F1=0
        for m in numpy.arange(1, 2 * self.ITERATIONS + 1, 2):
            F1 += 4. / (m * numpy.pi) * mpmath.sin(0.5 * m * numpy.pi * y_star) * (1. / (s + m**2 * numpy.pi**2))
        return F1

    def F2_LT(self, y_star, s):
        # m = numpy.arange(1, 2 * self.ITERATIONS + 1, 2)
        # F2 = numpy.zeros(y_star.shape)
        # np_cos = numpy.frompyfunc(mpmath.cos, 1, 1)
        # for i_y, y in enumerate(y_star):
        #     F2[i_y] = numpy.sum((8.0 / (m * numpy.pi) ** 2) * np_cos(0.5 * m * numpy.pi * y) * ((1.0 / s) - (1.0 / (s + (m * numpy.pi)**2))))
        # return F2
        F2 = 0.
        for m in numpy.arange(1, 2 * self.ITERATIONS + 1, 2):
            F2 += (8. / (m * numpy.pi)**2) * mpmath.cos(0.5 * m * numpy.pi * y_star) * (1. / s - (1. / (s + m**2 * numpy.pi**2)))
        return F2

    def F3_LT(self, y_star, s):
        # m = numpy.arange(1, 2 * self.ITERATIONS + 1, 2)
        # F3 = numpy.zeros(y_star.shape)
        # for i_y, y in enumerate(y_star):
        #     F3[i_y] = numpy.sum((-4.0 / (m * numpy.pi * H)) * mpmath.sin(0.5 * m * numpy.pi * y) * ((1.0 / s) - (1.0 / (s + (m * numpy.pi)**2))))
        # return F3

        F3 = 0.
        for m in numpy.arange(1, 2 * self.ITERATIONS + 1, 2):
            F3 += (-4.0 / (m * numpy.pi * H)) * mpmath.sin(0.5 * m * numpy.pi * y_star) * (1.0 / s - (1. / (s + m**2 * numpy.pi**2)))
        return F3

    def strain(self, locs):
        """Compute strain field at locations."""
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_zz = 0.0
        e_xy = 0.0

        strain = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[:, :, 0] = e_xx
        strain[:, :, 1] = e_yy
        strain[:, :, 2] = e_zz
        strain[:, :, 3] = e_xy
        return strain

    def stress(self, locs):
        """Compute stress field at locations."""
        (npts, dim) = locs.shape
        ntpts = tsteps.shape[0]
        poisson_ratio = (3 * K_d - 2 * G) / (2 * (3 * K_d + G))
        trace_strain = self.trace_strain(locs)
        pressure = self.pressure(locs)
        e_xx = 0.0
        e_yy = self.trace_strain(locs)
        e_xy = 0.0

        stress = numpy.zeros((ntpts, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[:, :, 0] = (
            ((2 * G * poisson_ratio) / (1 - 2 * poisson_ratio)) * trace_strain
            + 2 * G * e_xx
            - alpha * pressure
        )
        stress[:, :, 1] = (
            ((2 * G * poisson_ratio) / (1 - 2 * poisson_ratio)) * trace_strain
            + 2 * G * e_yy
            - alpha * pressure
        )
        stress[:, :, 2] = (
            (2 * G * poisson_ratio) / (1 - 2 * poisson_ratio)
        ) * trace_strain - alpha * pressure
        stress[:, :, 3] = 2 * G * e_xy
        return stress

    def bc_ypos_traction(self, locs):
        """Compute initial traction at locations."""
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:, :, 0] = 0.0
        traction[:, :, 1] = -P_0
        return traction

# End of file
