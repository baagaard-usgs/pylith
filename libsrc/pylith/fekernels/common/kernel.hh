// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#  define PYLITH_HOST_DEVICE __host__ __device__
#else
#  define PYLITH_HOST_DEVICE
#endif

// Shorthand for every kernel helper
#define PYLITH_KERNEL PYLITH_HOST_DEVICE inline
