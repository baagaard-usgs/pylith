# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from .InitializePhase import InitializePhase

from .topology import ReverseCuthillMcKee


class MeshReordering(InitializePhase):
    """
    Reorder mesh using the reverse Cuthill-McKee algorithm.

    Implements `InitializePhase`.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.mesh_initializer.reorder_mesh]
        """
    }

    def __init__(self, name="reorder_mesh"):
        """Constructor."""
        InitializePhase.__init__(self, name)

    def initialize(self, mesh, problem):
        """Do mesh initialize phase."""
        ReverseCuthillMcKee.reorder(mesh)
        return mesh

    def _configure(self):
        """Set members based on inventory."""
        InitializePhase._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////


def initialize_phase():
    """Factory associated with MeshReordering."""
    return MeshReordering()


# End of file
