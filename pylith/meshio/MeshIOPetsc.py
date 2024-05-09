# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import pathlib

from .MeshIOObj import MeshIOObj
from .meshio import MeshIOPetsc as ModuleMeshIOPetsc


class MeshIOPetsc(MeshIOObj, ModuleMeshIOPetsc):
    """
    Python object for a variety of reading/writing finite-element meshes using PETSc.
    Currently, the primary use of this object is to import meshes from Gmsh.

    :::{warning}
    The coordinate system associated with the mesh must be a Cartesian coordinate system, such as a generic Cartesian coordinate system or a geographic projection.
    :::

    Implements `MeshIOObj`.
    """

    import pythia.pyre.inventory

    filename = pythia.pyre.inventory.str("filename", default="")
    filename.meta['tip'] = "Name of mesh file for reading with PETSc."

    prefix = pythia.pyre.inventory.str("options_prefix", default="")
    prefix.meta['tip'] = "Name of PETSc options prefix for this mesh."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pythia.pyre.inventory.facility("coordsys", family="coordsys", factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."

    def __init__(self, mode=MeshIOObj.READ, name="meshiopetsc"):
        """Constructor.
        """
        MeshIOObj.__init__(self, mode, name)

    def preinitialize(self):
        """Do minimal initialization."""
        MeshIOObj.preinitialize(self)
        ModuleMeshIOPetsc.setFilename(self, self.filename)
        ModuleMeshIOPetsc.setPrefix(self, self.prefix)

    def _configure(self):
        """Set members based using inventory.
        """
        MeshIOObj._configure(self)

    def _createModuleObj(self):
        """Create C++ MeshIOPetsc object.
        """
        ModuleMeshIOPetsc.__init__(self)

    def _validate(self, context):
        if 0 == len(self.filename) and self.mode == self.READ and not pathlib.Path(self.filename).is_file():
            context.error(IOError(f"Input mesh '{self.filename}' not found."))

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_input():
    """Factory associated with MeshIOPetsc.
    """
    return MeshIOPetsc(mode=MeshIOPetsc.READ)


def mesh_output():
    """Factory associated with MeshIOPetsc.
    """
    return MeshIOPetsc(mode=MeshIOPetsc.WRITE)


# End of file
