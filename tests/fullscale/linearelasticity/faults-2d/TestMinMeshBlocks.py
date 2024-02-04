#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest

from pylith.testing.FullTestApp import (FullTestCase, Check)

import meshes
import minmesh_blocks_soln


# -------------------------------------------------------------------------------------------------
class TestCase(FullTestCase):

    def setUp(self):
        defaults = {
            "filename": "output/{name}-{mesh_entity}.h5",
            "exact_soln": minmesh_blocks_soln.AnalyticalSoln(),
            "mesh": self.mesh,
        }
        self.checks = [
            Check(
                mesh_entities=["domain", "boundary_ypos"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_xneg", "mat_xpos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                cell_fields = ["density", "bulk_modulus", "shear_modulus"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_xneg", "mat_xpos"],
                vertex_fields = ["displacement"],
                cell_fields = ["cauchy_strain"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["mat_xneg", "mat_xpos"],
                cell_fields = ["cauchy_stress"],
                scale = 1.0e+6,
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields=["initial_amplitude"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["bc_xneg", "bc_xpos"],
                vertex_fields=["displacement"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                filename="output/{name}-{mesh_entity}_info.h5",
                vertex_fields = ["normal_dir", "strike_dir"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                vertex_fields=["slip"],
                defaults=defaults,
            ),
            Check(
                mesh_entities=["fault"],
                vertex_fields=["traction_change"],
                defaults=defaults,
                scale = 1.0e+6
            ),
        ]

    def run_pylith(self, testName, args):
        FullTestCase.run_pylith(self, testName, args)


# -------------------------------------------------------------------------------------------------
class TestMinMeshQuad(TestCase):

    def setUp(self):
        self.name = "minmesh_blocks_quad"
        self.mesh = meshes.MinMeshQuad()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["pylithapp_minmesh.cfg", "minmesh_blocks.cfg", "minmesh_blocks_quad.cfg"])
        return


# -------------------------------------------------------------------------------------------------
class TestMinMeshTri(TestCase):

    def setUp(self):
        self.name = "minmesh_blocks_tri"
        self.mesh = meshes.MinMeshTri()
        super().setUp()

        TestCase.run_pylith(self, self.name, ["pylithapp_minmesh.cfg", "minmesh_blocks.cfg", "minmesh_blocks_tri.cfg"])
        return


# -------------------------------------------------------------------------------------------------
def test_cases():
    return [
        TestMinMeshQuad,
        TestMinMeshTri,
    ]


# -------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    FullTestCase.parse_args()

    suite = unittest.TestSuite()
    for test in test_cases():
        suite.addTest(unittest.makeSuite(test))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
