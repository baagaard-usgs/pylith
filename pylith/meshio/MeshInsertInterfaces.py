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


class MeshInsertInterfaces(InitializePhase):
    """
    Insert fault interfaces using PETSc transform operation.

    Implements `InitializePhase`.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.mesh_initializer.insert_interfaces]
        """
    }

    def __init__(self, name="mesh_insert_interfaces"):
        """Constructor."""
        InitializePhase.__init__(self, name)

    def initialize(self, mesh, problem):
        """Do mesh initialize phase."""
        cohesiveLabelValue = 100
        for material in problem.materials.components():
            labelValue = material.labelValue
            cohesiveLabelValue = max(cohesiveLabelValue, labelValue + 1)
        for interface in problem.interfaces:
            if mpi_is_root():
                self._info.log("Inserting fault interface '%s'." % interface.labelName)
            interface.preinitialize(problem)
            interface.setCohesiveLabelValue(cohesiveLabelValue)
            mesh = interface.transformTopology(mesh)
            cohesiveLabelValue += 1
        return mesh

    def _configure(self):
        """Set members based on inventory."""
        InitializePhase._configure(self)


# FACTORIES ////////////////////////////////////////////////////////////


def initialize_phase():
    """Factory associated with MeshInsertInterfaces."""
    return MeshInsertInterfaces()


# End of file
