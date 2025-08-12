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


class InitializePhase(Component):
    """
    Abstract base class for mesh initialization phases.
    """

    def __init__(self, name="initialize_phase"):
        """Constructor."""
        Component.__init__(self, name, facility="initialize_phase")

    def initialize(self, mesh, problem):
        """Do mesh initialization phase."""
        raise NotImplementedError()


# End of file
