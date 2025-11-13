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

#include "pylith/problems/problemsfwd.hh" // forward declarations
#include "pylith/topology/SubfieldFactory.hh" // ISA SubfieldFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES TimeHistory

class pylith::problems::SolutionFactory : public pylith::topology::SubfieldFactory {
    friend class TestSolutionFactory; // unit testing

    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    static const std::string displacement;
    static const std::string velocity;
    static const std::string pressure;
    static const std::string fluid_pressure;
    static const std::string fluid_pressure_dot;
    static const std::string trace_strain;
    static const std::string trace_strain_dot;
    static const std::string lagrange_multiplier_fault;
    static const std::string temperature;

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    SolutionFactory(void);

    /// Destructor.
    ~SolutionFactory(void);

    /** Set values from spatial database.
     *
     * @param[in] db Spatial database for solution values.
     */
    void setValues(const std::shared_ptr<spatialdata::spatialdb::SpatialDB>& db);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get subfield description.
     *
     * @param[in] subfieldName Name of subfield.
     * @returns Description of subfield.
     */
    pylith::topology::FieldBase::Description _getDescription(const std::string& subfieldName) const override;

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    SolutionFactory(const SolutionFactory &) = delete;
    const SolutionFactory& operator=(const SolutionFactory&) = delete;

}; // class SolutionFactory

// End of file
