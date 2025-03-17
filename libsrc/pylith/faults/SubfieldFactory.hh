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

#include "pylith/faults/faultsfwd.hh" // forward declarations
#include "pylith/topology/SubfieldFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::SubfieldFactory : public pylith::topology::SubfieldFactory {
    friend class TestAuxiliaryFieldFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    SubfieldFactory(void);

    /// Destructor.
    ~SubfieldFactory(void);

    /// Add slip subfield to auxiliary field.
    void addSlip(void);

    /// Add slip rate subfield to auxiliary field.
    void addSlipRate(void);

    /// Add slip acceleration subfield to auxiliary field.
    void addSlipAcceleration(void);

    /// Add traction change subfield to derived field.
    void addTractionChange(void);

    /// Add fault normal direction subfield to auxiliary field.
    void addNormalDir(void);

    /// Add fault strike direction subfield to auxiliary field.
    void addStrikeDir(void);

    /// Add fault up-dip direction modulus subfield to auxiliary field.
    void addUpDipDir(void);

    /// Add subfields using discretizations provided.
    virtual void addSubfields(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    SubfieldFactory(const SubfieldFactory &); ///< Not implemented.
    const SubfieldFactory& operator=(const SubfieldFactory&); ///< Not implemented

}; // class SubfieldFactory

// End of file
