#!/usr/bin/env nemesis

# Import Gmsh Python interface
import gmsh
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh)

class App(GenerateMesh):
    """
    Domain is a 30km x 30 km square.
    Outer shell is a 5 km radius circle.
    Innter shell is a 1 km radius circle.
    """
    DOMAIN_X = DOMAIN_Y = 30e3
    INNER_SHELL_RAD = 1e3
    OUTER_SHELL_RAD = 5e3
    SHELL_CENTER = (DOMAIN_X/2, -DOMAIN_Y/2, 0)
    DX_CHAMBER = 250
    DX_BIAS = 1.1

    def __init__(self):
        # Set the cell choices available through command line options.
        # The default cell type `tri` and filename match the mesh used
        # in the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "pve_mesh_tri.msh"

    def create_geometry(self):
        """Create geometry.
        """
        # Create domain surface
        domain_box = gmsh.model.occ.add_rectangle(0, -self.DOMAIN_Y, 0, self.DOMAIN_X, self.DOMAIN_Y)
        gmsh.model.occ.add_plane_surface([domain_box])

        # Create inner chamber surface
        inner_chamber_ellipse = gmsh.model.occ.add_ellipse(self.SHELL_CENTER[0], self.SHELL_CENTER[1], self.SHELL_CENTER[2], self.INNER_SHELL_RAD, self.INNER_SHELL_RAD)
        inner_chamber_loop = gmsh.model.occ.add_curve_loop([inner_chamber_ellipse])
        self.inner_chamber_surface = gmsh.model.occ.add_plane_surface([inner_chamber_loop])

        # Create outer chamber surface
        outer_chamber_ellipse = gmsh.model.occ.add_ellipse(self.SHELL_CENTER[0], self.SHELL_CENTER[1], self.SHELL_CENTER[2], self.OUTER_SHELL_RAD, self.OUTER_SHELL_RAD)
        outer_chamber_loop = gmsh.model.occ.add_curve_loop([outer_chamber_ellipse])
        outer_chamber_surface_temp = gmsh.model.occ.add_plane_surface([outer_chamber_loop])
        
        # embed inner chamber to outer chamber
        self.outer_chamber_surface = gmsh.model.occ.fragment([(2, outer_chamber_surface_temp)], [(2, self.inner_chamber_surface)])[0][1][1]

        # embed chambers to domain
        self.domain_surface = gmsh.model.occ.fragment([(2, domain_box)], [(2, self.outer_chamber_surface)])[0][1][1]

        # get rid of overlap
        gmsh.model.occ.remove_all_duplicates()

        # Get boundaries
        _, domain_curves = gmsh.model.occ.get_curve_loops(self.domain_surface)
        self.y_neg_boundary = domain_curves[0][0]
        self.x_neg_boundary = domain_curves[0][3]
        self.x_pos_boundary = domain_curves[0][1]
        self.y_pos_boundary = domain_curves[0][2]

        _, outer_chamber_curves = gmsh.model.occ.get_curve_loops(self.outer_chamber_surface)
        self.outer_chamber_boundary = outer_chamber_curves[0][0]

        _, inner_chamber_curves = gmsh.model.occ.get_curve_loops(self.inner_chamber_surface)
        self.inner_chamber_boundary = inner_chamber_curves[0][0]

        # dim_tags = gmsh.model.occ.get_entities(dim=2)
        # self.boundaries = [tag for dim, tag in dim_tags]

        gmsh.model.occ.synchronize()

    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.
        """

        # Create materials for the domain, outer chamber and inner chamber. 
        # The tag argument specifies the integer tag for the physical group.
        # The entities argument specifies the array of surfaces for that material.
        materials = (
            MaterialGroup(tag=1, entities=[self.domain_surface]),
            MaterialGroup(tag=2, entities=[self.outer_chamber_surface]),
            MaterialGroup(tag=3, entities=[self.inner_chamber_surface])
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries.
        # We use the `BoundaryGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.
        boundary_groups = (
            VertexGroup(name="boundary_xneg", tag=1, dim=1, entities=[self.x_neg_boundary]),
            VertexGroup(name="boundary_xpos", tag=2, dim=1, entities=[self.x_pos_boundary]),
            VertexGroup(name="boundary_yneg", tag=3, dim=1, entities=[self.y_neg_boundary]),
            VertexGroup(name="boundary_ypos", tag=4, dim=1, entities=[self.y_pos_boundary]),
            VertexGroup(name="boundary_outer_chamber", tag=5, dim=1, entities=[self.outer_chamber_boundary]),
            VertexGroup(name="boundary_inner_chamber", tag=6, dim=1, entities=[self.inner_chamber_boundary]),
        )
        for group in boundary_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh. Should also include optimizing the mesh quality.
        """
        # Set discretization size with geometric progression from distance to the chamber boundary.

        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the reservoir boundary.
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "CurvesList", [self.inner_chamber_boundary, self.outer_chamber_boundary])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the boundary, the distance from
        # the boundary (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression`
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(distance, min_dx=self.DX_CHAMBER, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        # Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        if cell == "quad":
            # Generate a tri mesh and then recombine cells to form quadrilaterals.
            # We use the Frontal-Delaunay for Quads algorithm.
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.recombine()
            gmsh.model.mesh.generate(2)
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")

if __name__ == "__main__":
    App().main()


# End of file
