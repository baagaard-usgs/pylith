# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================

from pythia.pyre.components.Component import Component


class SerialPhases(Component):
    """
    Mesh initialization phases for reading in serial.

    :::{seealso}
    See [`MeshInitializer` Component](../problems/Problem.md).
    :::
    """

    import pythia.pyre.inventory

    from .MeshReader import MeshReader

    read_mesh = pythia.pyre.inventory.facility(
        "read_mesh", family="initialize_phase", factory=MeshReader
    )
    read_mesh.meta["tip"] = "Read mesh in serial."

    from .MeshReordering import MeshReordering

    reorder_mesh = pythia.pyre.inventory.facility(
        "reorder_mesh", family="initialize_phase", factory=MeshReordering
    )
    reorder_mesh.meta["tip"] = "Reorder mesh using reverse Cuthill-McKee algorithm."

    from .MeshDistributor import MeshDistributor

    distribute_mesh = pythia.pyre.inventory.facility(
        "distribute_mesh", family="initialize_phase", factory=MeshDistributor
    )
    distribute_mesh.meta["tip"] = "Distribute mesh among processes."

    from .MeshInsertInterfaces import MeshInsertInterfaces

    insert_interfaces = pythia.pyre.inventory.facility(
        "insert_interfaces", family="initialize_phase", factory=MeshInsertInterfaces
    )
    insert_interfaces.meta["tip"] = (
        "Insert interfaces using PETSc mesh transform operation."
    )

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="serial_phases"):
        """Constructor."""
        Component.__init__(self, name, facility="serial_phases")

    def components(self):
        """Order of facilities in Inventory is ambiguous, so overwrite
        components() to insure order is [displacement, pressure, trace_strain].

        """
        return [
            self.read_mesh,
            self.reorder_mesh,
            self.distribute_mesh,
            self.insert_interfaces,
        ]


def phasesFactory(name):
    """Factory for output items."""
    from pythia.pyre.inventory import facility
    from pylith.meshio.InitializePhase import InitializePhase

    return facility(name, family="initialize_phase", factory=InitializePhase)


class MeshInitializer(Component):
    """
    Manager for reading and setting up a finite-element mesh.
    """

    DOC_CONFIG = {
        "cfg": """
            [pylithapp.mesh_initializer]
            phases = [read_mesh, reorder_mesh, distribute_mesh, insert_interfaces]
            read_mesh = pylith.meshio.MeshReader
            reorder_mesh = pylith.meshio.MeshReordering
            distribute_mesh = pylith.meshio.MeshDistributor
            insert_interfaces = pylith.meshio.MeshInsertInterfaces
        """
    }

    import pythia.pyre.inventory

    phases = pythia.pyre.inventory.facilityArray(
        "phases", itemFactory=phasesFactory, factory=SerialPhases
    )
    phases.meta["tip"] = "Phases for mesh initialization."

    def __init__(self, name="mesh_initializer"):
        """Constructor."""
        Component.__init__(self, name, facility="mesh_initializer")

    def initialize(self):
        """Read and setup a finite-element mesh."""
        mesh = None
        for phase in self.phases:
            mesh = phase.initializer(mesh)
        return mesh

    def _configure(self):
        """Set members based using inventory."""
        Component._configure(self)


# End of file
