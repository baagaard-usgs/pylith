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

MMS_DRIVER = mmstest_linearporoelasticity_nofaults2d

TESTS = \
	run_PressureGradient.sh \
	run_PressureGradientPVE.sh \
	run_TerzaghiPVE.sh \
	run_SimpleShearPVE.sh

check_SCRIPTS = $(TESTS)

check_PROGRAMS = $(MMS_DRIVER)

mmstest_linearporoelasticity_nofaults2d_SOURCES = \
	$(top_srcdir)/tests/src/MMSTest.cc \
	$(top_srcdir)/tests/src/driver_catch2.cc \
	TestLinearPoroelasticity.cc \
	TestLinearPoroviscoelasticity.cc \
	TestCases.cc \
	PressureGradient.cc \
	PressureGradientPVE.cc \
	TerzaghiPVE.cc \
	SimpleShearPVE.cc

run_%.sh:
	echo "#!/bin/bash" > $@
	echo "$(abs_builddir)/$(MMS_DRIVER) [$*]" >> $@
	chmod +x $@

dist_noinst_HEADERS = \
	TestLinearPoroelasticity.hh	\
	TestLinearPoroviscoelasticity.hh \
	PressureGradient.hh \
	PressureGradientPVE.hh \
	TerzaghiPVE.hh \
	SimpleShearPVE.hh


dist_noinst_DATA = \
	data/tri.mesh \
	data/quad.mesh \
	data/quad_simple.mesh

export_datadir = $(abs_builddir)
include $(top_srcdir)/tests/data.am


# End of file
