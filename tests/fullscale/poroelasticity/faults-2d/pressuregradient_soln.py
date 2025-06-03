# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @brief Analytical solution to pressure gradient for poroelasticity with zero prescribed slip.
#

import numpy


# Physical properties
p_density = 2500.0

p_mu = p_density * 0**2
p_lambda = p_density * 0**2 - 2 * p_mu


# ----------------------------------------------------------------------
class AnalyticalSoln(object):
    """Analytical solution to shear problem.
    """
    SPACE_DIM = 2
    TENSOR_SIZE = 4

    def __init__(self):
        self.fields = {
            "displacement": self.displacement,
            "density": self.density,
            "shear_modulus": self.shear_modulus,
            "solid_bulk_modulus": self.solid_bulk_modulus,
            "cauchy_strain": self.strain,
            "cauchy_stress": self.stress,
            "initial_amplitude": {
                "bc_disp_xneg": self.displacement,
                "bc_disp_xpos": self.displacement,
                "bc_disp_yneg": self.displacement,
                "bc_disp_ypos": self.displacement,
                "bc_press_xneg": self.displacement,
                "bc_press_xpos": self.displacement,
            },
            "normal_dir": {
                "bc_disp_xneg": self.orientation_dir((-1, 0)),
                "bc_disp_xpos": self.orientation_dir((+1, 0)),
                "bc_disp_yneg": self.orientation_dir((0, -1)),
                "bc_disp_ypos": self.orientation_dir((0, +1)),
                "bc_press_xneg": self.orientation_dir((-1, 0)),
                "bc_press_xpos": self.orientation_dir((+1, 0)),
                "fault": self.orientation_dir((+1, 0)),
            },
            "tangential_dir": {
                "bc_disp_xneg": self.orientation_dir((0, -1)),
                "bc_disp_xpos": self.orientation_dir((0, +1)),
                "bc_disp_yneg": self.orientation_dir((+1, 0)),
                "bc_disp_ypos": self.orientation_dir((-1, 0)),
                "bc_press_xneg": self.orientation_dir((0, -1)),
                "bc_press_xpos": self.orientation_dir((0, +1)),
            },
            "slip": self.slip,
            "traction_change": self.traction_change,
            "strike_dir": self.orientation_dir((0, +1)),
        }

    def getField(self, name, mesh_entity, pts):
        if isinstance(self.fields[name], dict):
            field = self.fields[name][mesh_entity](pts)
        else:
            field = self.fields[name](pts)
        return field

    def displacement(self, locs):
        """Compute displacement field at locations.
        """
        (npts, dim) = locs.shape
        disp = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        disp[0,:, 0] = 0 * locs[:, 0] + 0 * locs[:, 1]
        disp[0,:, 1] = 0 * locs[:, 1] + 0 * locs[:, 0]
        return disp

    def density(self, locs):
        """Compute density field at locations.
        """
        (npts, dim) = locs.shape
        density = p_density * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return density

    def shear_modulus(self, locs):
        """Compute shear modulus field at locations.
        """
        (npts, dim) = locs.shape
        shear_modulus = p_mu * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return shear_modulus

    def solid_bulk_modulus(self, locs):
        """Compute bulk modulus field at locations.
        """
        (npts, dim) = locs.shape
        bulk_modulus = (p_lambda + 2.0 / 3.0 * p_mu) * numpy.ones((1, npts, 1), dtype=numpy.float64)
        return bulk_modulus

    def strain(self, locs):
        """Compute strain field at locations.
        """
        (npts, dim) = locs.shape
        strain = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        strain[0,:, 0] = 0
        strain[0,:, 1] = 0
        strain[0,:, 2] = 0
        strain[0,:, 3] = 0
        return strain

    def stress(self, locs):
        """Compute stress field at locations.
        """
        (npts, dim) = locs.shape
        stress = numpy.zeros((1, npts, self.TENSOR_SIZE), dtype=numpy.float64)
        stress[0,:, 0] = 0
        stress[0,:, 1] = 0
        stress[0,:, 2] = 0
        stress[0,:, 3] = 0
        return stress

    def slip(self, locs):
        """Compute slip field.
        """
        (npts, dim) = locs.shape
        slip = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        return slip

    def traction_change(self, locs):
        """Compute change in traction on faults.
        """
        (npts, dim) = locs.shape
        traction = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
        traction[:,:,1] = 0
        return traction

    def orientation_dir(self, vector):
        def fn_dir(locs):
            (npts, dim) = locs.shape
            values = numpy.zeros((1, npts, self.SPACE_DIM), dtype=numpy.float64)
            for d in range(self.SPACE_DIM):
                values[:,:,d] = vector[d]
            return values
        return fn_dir


# End of file
