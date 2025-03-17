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

#include "pylith/bc/bcfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include <string> // HASA std::string

class pylith::bc::BoundaryCondition : public pylith::problems::Physics {
    friend class TestBoundaryCondition; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    BoundaryCondition(void);

    /// Destructor.
    virtual ~BoundaryCondition(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set name of solution subfield associated with boundary condition.
     *
     * @param[in] value Name of solution subfield.
     */
    void setSubfieldName(const char* value);

    /** Get name of solution subfield associated with boundary condition.
     *
     * @preturn Name of solution subfield.
     */
    const char* getSubfieldName(void) const;

    /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void setRefDir1(const pylith::real vec[3]);

    /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void setRefDir2(const pylith::real vec[3]);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const override;

    /** Create diagnostic field.
     *
     * @param[in] solution Solution field.
     * @param[in] physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Diagnostic field if applicable, otherwise NULL.
     */
    virtual
    pylith::topology::Field* createDiagnosticField(const pylith::topology::Field& solution,
                                                   const pylith::topology::Mesh& physicsMesh) override;

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get subfield factory associated with physics.
     *
     * @return Subfield factory for physics object.
     */
    pylith::topology::SubfieldFactory* _getSubfieldFactory(void) override;

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const pylith::real dt);

    /** Set kernels for computing diagnostic field.
     *
     * @param[out] integrator Integrator for material.
     * @param[in] solution Solution field.
     */
    virtual
    void _setKernelsDiagnosticField(pylith::feassemble::IntegratorBoundary* integrator,
                                    const pylith::topology::Field& solution) const;

    /** Set kernels for computing diagnostic field.
     *
     * @param[out] constraint Integrator for material.
     * @param[in] solution Solution field.
     */
    virtual
    void _setKernelsDiagnosticField(pylith::feassemble::Constraint* constraint,
                                    const pylith::topology::Field& solution) const;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::real _refDir1[3]; ///< First choice reference direction used to compute boundary tangential directions.
    pylith::real _refDir2[3]; ///< Second choice reference direction used to compute boundary tangential directions.
    std::string _subfieldName; ///< Name of solution subfield for boundary condition.
    std::unique_ptr<pylith::bc::SubfieldFactory> _subfieldFactory; ///< Factory for auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    BoundaryCondition(const BoundaryCondition&); ///< Not implemented.
    const BoundaryCondition& operator=(const BoundaryCondition&); ///< Not implemented.

}; // class BoundaryCondition

// End of file
