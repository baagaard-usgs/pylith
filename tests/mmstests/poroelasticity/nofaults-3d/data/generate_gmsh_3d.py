#!/usr/bin/env nemesis

import gmsh
from pylith.meshio.gmsh_utils import (BoundaryGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Block is DOMAIN_X by DOMAIN_Y by DOMAIN_Z with discretization size DX.
    """
    DOMAIN_X = DOMAIN_Y = DOMAIN_Z = 8.0e+3
    DX = 2.0e+3

    def __init__(self):
        self.cell_choices = {
            "required": True,
            "choices": ["tet", "hex"],
            }
        self.filename = "tet.msh"

    def create_geometry(self):
        """Create geometry.
        """
        self.v_domain = gmsh.model.occ.add_box(-0.5*self.DOMAIN_X, -0.5*self.DOMAIN_Y, -self.DOMAIN_Z,
                                               self.DOMAIN_X, self.DOMAIN_Y, self.DOMAIN_Z)

        gmsh.model.occ.synchronize()

        dim_tags = gmsh.model.occ.get_entities(dim=2)
        self.boundaries = [tag for dim, tag in dim_tags]

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """
        materials = (
            MaterialGroup(tag=24, entities=[self.v_domain]),
        )
        for material in materials:
            material.create_physical_group()

        boundary_groups = (
            BoundaryGroup(name="boundary_xneg", tag=1, dim=2, entities=[1]),
            BoundaryGroup(name="boundary_xpos", tag=2, dim=2, entities=[2]),
            BoundaryGroup(name="boundary_yneg", tag=3, dim=2, entities=[3]),
            BoundaryGroup(name="boundary_ypos", tag=4, dim=2, entities=[4]),
            BoundaryGroup(name="boundary_zpos", tag=5, dim=2, entities=[6]),
            BoundaryGroup(name="boundary_zneg", tag=6, dim=2, entities=[5]),
        )
        for group in boundary_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh. Should also include optimizing the mesh quality.
        """
        gmsh.option.setNumber("Mesh.MeshSizeMin", self.DX)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.DX)
        if cell == "hex":
            gmsh.model.mesh.set_transfinite_automatic(recombine=True)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
