# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

include $(top_srcdir)/tests/check_catch2.am

MMS_DRIVER = mmstest_linearporoelasticity_nofaults3d

TESTS = \
	run_SimpleShearPVE3D.sh \
	run_TerzaghiPVE3D.sh

check_SCRIPTS = $(TESTS)

check_PROGRAMS = $(MMS_DRIVER)

mmstest_linearporoelasticity_nofaults3d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestLinearPoroviscoelasticity.cc \
	TestCases.cc \
	SimpleShearPVE3D.cc \
	TerzaghiPVE3D.cc

run_%.sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
	chmod +x $@

dist_noinst_HEADERS = \
	TestLinearPoroviscoelasticity.hh \
	SimpleShearPVE3D.hh \
	TerzaghiPVE3D.hh


dist_noinst_DATA = \
	data/hex.msh \
	data/tet.msh \
	data/tet_all_boundaries.msh \
	data/hex_all_boundaries.msh

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
