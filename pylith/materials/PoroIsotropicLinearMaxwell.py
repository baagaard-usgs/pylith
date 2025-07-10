# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pylith/materials/IsotropicLinearPoroelasticity.py
#
# @brief Python material for isotropic, linearly elastic, plane
# strain material.
#
# Factory: poroelasticity_rheology

from .RheologyPoroelasticity import RheologyPoroelasticity
from .materials import PoroIsotropicLinearMaxwell as ModulePoroIsotropicLinearMaxwell


class PoroIsotropicLinearMaxwell(RheologyPoroelasticity, ModulePoroIsotropicLinearMaxwell):
    """
    Poroelasticity bulk rheology.

    Implements `RheologyPoroelasticity`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.materials.mat_incompelastic.rheology]
            use_reference_state = False

            auxiliary_subfields.shear_modulus.basis_order = 0
            auxiliary_subfields.bulk_modulus.basis_order = 0
        """
    }

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    useTensorPermeability = pythia.pyre.inventory.bool("use_tensor_permeability", default=False)
    useTensorPermeability.meta['tip'] = "Use tensor permeability."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="poroisotropiclinearmaxwell"):
        """Constructor.
        """
        RheologyPoroelasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsPoroIsotropicLinearMaxwell import AuxSubfieldsPoroIsotropicLinearMaxwell
        self.auxiliarySubfields = AuxSubfieldsPoroIsotropicLinearMaxwell("auxiliary_subfields")

        from .DerivedSubfieldsPoroelasticity import DerivedSubfieldsPoroelasticity
        self.derivedSubfields = DerivedSubfieldsPoroelasticity("derived_subfields")

    def preinitialize(self, mesh):
        RheologyPoroelasticity.preinitialize(self, mesh)

        ModulePoroIsotropicLinearMaxwell.useReferenceState(self, self.useReferenceState)
        ModulePoroIsotropicLinearMaxwell.useTensorPermeability(self, self.useTensorPermeability)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModulePoroIsotropicLinearMaxwell.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def poroviscoelasticity_rheology():
    """Factory associated with PoroIsotropicLinearMaxwell.
    """
    return PoroIsotropicLinearMaxwell()


# End of file
