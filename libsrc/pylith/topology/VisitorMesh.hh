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

#include "pylith/topology/topologyfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h" // HASA PetscVec, PetscSection
#include "pylith/utils/arrayfwd.hh" // USES scalar_array

// VecVisitorMesh ----------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
class pylith::topology::VecVisitorMesh { // VecVisitorMesh
    friend class TestVecVisitorMesh; // unit testing

    // PUBLIC ENUMS ///////////////////////////////////////////////////////////////////////////////
public:

    enum SectionEnum {
        LOCAL_SECTION=0,
        GLOBAL_SECTION=1,
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor with field over a mesh.
     *
     * The optional subfield argument is designed to improve performance
     * when the visitor is associated with a single subfield within a
     * field.
     *
     * @param field Field over a mesh.
     * @param subfield Name of subfield section to use instead of field section.
     */
    VecVisitorMesh(const Field& field,
                   const char* subfield=0);

    /// Default destructor
    ~VecVisitorMesh(void);

    /** Initialize using field over a mesh or submesh.
     *
     * The optional subfield argument is designed to improve performance
     * when the visitor is associated with a single subfield within a
     * field.
     *
     * @param field Field over a mesh/submesh.
     * @param subfield Name of subfield section to use instead of field section.
     */
    void initialize(const Field& field,
                    const char *subfield=0);

    /// Clear cached data.
    void clear(void);

    /** Set selection (local or global) of section.
     *
     * @param[in] value
     */
    void selectSection(pylith::topology::VecVisitorMesh::SectionEnum value);

    /** Get the PETSc current section.
     *
     * @returns PETSc current section.
     */
    PetscSection selectedSection(void) const;

    /** Get the array of values associated with the local PETSc Vec.
     *
     * @returns Array of values.
     */
    pylith::scalar* localArray(void) const;

    /** Get the local PETSc Vec.
     *
     * @returns PETSc Vec.
     */
    PetscVec localVec(void) const;

    /** Get fiber dimension for values at point.
     *
     * @param point Point in mesh.
     * @returns Fiber dimension.
     */
    pylith::integer sectionDof(const pylith::integer point) const;

    /** Get fiber dimension for values of subfield at point.
     *
     * @param subfield Index of subfield in section.
     * @param point Point in mesh.
     * @returns Fiber dimension.
     */
    pylith::integer sectionSubfieldDof(const pylith::integer subfield,
                                       const pylith::integer point) const;

    /** Get offset into values array for point.
     *
     * @param point Point in mesh.
     * @returns Offset.
     */
    pylith::integer sectionOffset(const pylith::integer point) const;

    /** Get offset into values array for point.
     *
     * @param subfield Index of subfield in section.
     * @param point Point in mesh.
     * @returns Offset.
     */
    pylith::integer sectionSubfieldOffset(const pylith::integer subfield,
                                          const pylith::integer point) const;

    /** Get array of values associated with closure.
     *
     * @pre Must be followed by call to restoreClosure().
     *
     * @param valuesCell Array of values for cell.
     * @param valuesSize Size of values array.
     * @param cell Finite-element cell.
     */
    void getClosure(pylith::scalar** valuesCell,
                    pylith::integer* valuesSize,
                    const pylith::integer cell) const;

    /** Get array of values associated with closure.
     *
     * @param values Array of values for cell.
     * @param cell Finite-element cell.
     */
    void getClosure(scalar_array* values,
                    const pylith::integer cell) const;

    /** Restore array of values associated with closure.
     *
     * @pre Must be preceded by call to getClosure().
     *
     * @param valuesCell Array of values for cell.
     * @param valuesSize Size of values array.
     * @param cell Finite-element cell.
     */
    void restoreClosure(pylith::scalar** valuesCell,
                        pylith::integer* valuesSize,
                        const pylith::integer cell) const;

    /** Set values associated with closure.
     *
     * @param valuesCell Array of values for cell.
     * @param valuesSize Size of values array.
     * @param cell Finite-element cell.
     * @param mode Mode for inserting values.
     */
    void setClosure(const pylith::scalar* valuesCell,
                    const pylith::integer valuesSize,
                    const pylith::integer cell,
                    const InsertMode mode) const;

    /** Optimize the closure operator by creating index for closures.
     *
     * :TODO: Remove this method. Call static version when setting up fields.
     */
    void optimizeClosure(void);

    /** Optimize the closure operator by creating index for closures.
     *
     * @param field Field to optimize closure for.
     */
    static
    void optimizeClosure(const Field& field);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    PetscDM _dm; ///< Cached PETSc dm for mesh.
    PetscVec _localVec; ///< Cached local PETSc Vec.
    PetscSection _localSection; ///< Cached PETSc local section.
    PetscSection _globalSection; ///< Cached PETSc global section.
    pylith::scalar* _localArray; ///< Cached local array
    SectionEnum _selectedSection; ///< Current selected section.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    VecVisitorMesh(const VecVisitorMesh&) = delete;
    const VecVisitorMesh& operator=(const VecVisitorMesh&) = delete;

};

// VecVisitorMesh

// MatVisitorMesh ----------------------------------------------------------
/** @brief Helper class for accessing field values at points in a
 *  finite-element mesh.
 */
class pylith::topology::MatVisitorMesh { // MatVisitorMesh
    friend class TestMatVisitorMesh; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param mat PETSc matrix.
     * @param field Field associated with matrix layout.
     */
    MatVisitorMesh(const PetscMat mat,
                   const Field& field);

    /// Default destructor
    ~MatVisitorMesh(void);

    // Initialize.
    void initialize(void);

    /// Clear cached data.
    void clear(void);

    /** Set values associated with closure.
     *
     * @param valuesCell Array of values for cell.
     * @param valuesSize Size of values array.
     * @param cell Finite-element cell.
     * @param mode Mode for inserting values.
     */
    void setClosure(const pylith::scalar* valuesCell,
                    const pylith::integer valuesSize,
                    const pylith::integer cell,
                    const InsertMode mode) const;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    const PetscMat _mat; ///< Cached PETSc matrix.
    PetscDM _dm; ///< Cached PETSc dm for mesh.
    PetscSection _localSection; ///< Cached PETSc local section.
    PetscSection _globalSection; ///< Cached PETSc global section.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    MatVisitorMesh(const MatVisitorMesh&) = delete;
    const MatVisitorMesh& operator=(const MatVisitorMesh&) = delete;

};

// MatVisitorMesh

#include "VisitorMesh.icc"

// End of file
