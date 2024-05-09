// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/meshio/MeshIO.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Stratum.hh" // USES Stratum

#include "pylith/utils/array.hh" // USES scalar_array, int_array
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_INFO
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept>

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIO::MeshIO(void) :
    _mesh(NULL) {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIO::~MeshIO(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIO::deallocate(void) {
} // deallocate


// ----------------------------------------------------------------------
// Get spatial dimension of mesh.
int
pylith::meshio::MeshIO::getMeshDim(void) const {
    return (_mesh) ? _mesh->getDimension() : 0;
} // getMeshDim


// ----------------------------------------------------------------------
// Read mesh from file.
void
pylith::meshio::MeshIO::read(pylith::topology::Mesh* mesh,
                             const bool checkTopology) {
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(!_mesh);

    _mesh = mesh;
    _read();

    PetscErrorCode err = PETSC_SUCCESS;

    // Check for bounding box with positive volume.
    PylithReal cmin[3];
    PylithReal cmax[3];
    err = DMGetBoundingBox(_mesh->getDM(), cmin, cmax);
    const PetscInt dim = _mesh->getDimension();
    PylithReal volume = 1.0;
    for (int i = 0; i < dim; ++i) {
        volume *= cmax[i] - cmin[i];
    } // for
    std::ostringstream msg;
    msg << "Domain bounding box:";
    for (int i = 0; i < dim; ++i) {
        msg << "\n    (" << cmin[i] << ", " << cmax[i] << ")";
    } // for
    PYLITH_COMPONENT_INFO_ROOT(msg.str());
    const PetscReal tolerance = 1.0e-8;
    if (volume < tolerance) {
        msg.clear();
        msg << "Domain bounding box volume (" << volume << ") is less than minimum tolerance ("
            << tolerance << "). This usually means you are trying to use a 2D mesh in 3D. Check that you are exporting "
            << "your mesh from the mesh generation software correctly and that your have specified the correct "
            << " coordinate system for the problem.";
        throw std::runtime_error(msg.str());
    } // if

    // Check mesh consistency
    if (checkTopology) {
        pylith::topology::MeshOps::checkTopology(*_mesh);
    } // if

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        _mesh->view("::ascii_info_detail");
    } // if
    // Respond to PETSc diagnostic output
    err = DMViewFromOptions(_mesh->getDM(), NULL, "-pylith_dm_view");PYLITH_CHECK_ERROR(err);

    _mesh = NULL;

    PYLITH_METHOD_END;
} // read


// ----------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIO::write(pylith::topology::Mesh* const mesh) { // write
    PYLITH_METHOD_BEGIN;

    assert(mesh);
    assert(!_mesh);

    _mesh = mesh;
    _write();
    _mesh = 0;

    PYLITH_METHOD_END;
} // write


// ----------------------------------------------------------------------
// Get coordinates of vertices in mesh.
void
pylith::meshio::MeshIO::_getVertices(scalar_array* coordinates,
                                     int* numVertices,
                                     int* spaceDim) const {
    PYLITH_METHOD_BEGIN;

    assert(coordinates);
    assert(numVertices);
    assert(spaceDim);
    assert(_mesh);

    *spaceDim = _mesh->getDimension();

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    PetscVec coordVec = NULL;
    PetscScalar* coordArray = NULL;
    PetscInt coordSize = 0;
    PylithScalar lengthScale = 1.0;
    PetscErrorCode err = 0;

    // Get length scale for dimensioning
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);

    // Get coordinates and dimensionalize values
    err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);
    err = VecGetArray(coordVec, &coordArray);PYLITH_CHECK_ERROR(err);
    err = VecGetLocalSize(coordVec, &coordSize);PYLITH_CHECK_ERROR(err);
    assert(coordSize % *spaceDim == 0);
    *numVertices = coordSize / *spaceDim;

    coordinates->resize(coordSize);
    for (PetscInt i = 0; i < coordSize; ++i) {
        (*coordinates)[i] = coordArray[i]*lengthScale;
    } // for
    err = VecRestoreArray(coordVec, &coordArray);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _getVertices


// ----------------------------------------------------------------------
// Get cells in mesh.
void
pylith::meshio::MeshIO::_getCells(int_array* cells,
                                  int* numCells,
                                  int* numCorners,
                                  int* meshDim) const {
    PYLITH_METHOD_BEGIN;

    assert(cells);
    assert(numCells);
    assert(numCorners);
    assert(meshDim);
    assert(_mesh);

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    *numCells = pylith::topology::MeshOps::getNumCells(*_mesh);assert(*numCells > 0);
    *numCorners = pylith::topology::MeshOps::getNumCorners(*_mesh);assert(*numCorners > 0);
    *meshDim = _mesh->getDimension();
    assert(cellsStratum.size() == *numCells);

    cells->resize((*numCells)*(*numCorners));

    PetscIS globalVertexNumbers = NULL;
    const PetscInt* gvertex = NULL;
    PetscErrorCode err = 0;

    err = DMPlexGetVertexNumbering(dmMesh, &globalVertexNumbers);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(globalVertexNumbers, &gvertex);PYLITH_CHECK_ERROR(err);
    for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
        DMPolytopeType ct;
        PetscInt nC = 0, closureSize, *closure = NULL;

        err = DMPlexGetCellType(dmMesh, c, &ct);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (PetscInt cl = 0; cl < closureSize*2; cl += 2) {
            if ((closure[cl] >= vStart) && (closure[cl] < vEnd)) {
                const PetscInt gv = gvertex[closure[cl]-vStart];
                (*cells)[index++] = gv < 0 ? -(gv+1) : gv;
                ++nC;
            }
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, c, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        err = DMPlexInvertCell(ct, &(*cells)[index-nC]);PYLITH_CHECK_ERROR(err);
        assert(nC == *numCorners);
    } // for

    PYLITH_METHOD_END;
} // _getCells


// ----------------------------------------------------------------------
// Tag cells in mesh with material identifiers.
void
pylith::meshio::MeshIO::_setMaterials(const int_array& materialIds) {
    PYLITH_METHOD_BEGIN;
    assert(_mesh);

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    PetscErrorCode err = 0;
    const char* const labelName = pylith::topology::Mesh::cells_label_name;

    if (!_mesh->getCommRank()) {
        PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
        topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
        const PetscInt cStart = cellsStratum.begin();
        const PetscInt cEnd = cellsStratum.end();

        if (size_t(cellsStratum.size()) != materialIds.size()) {
            std::ostringstream msg;
            msg << "Mismatch in size of materials identifier array ("
                << materialIds.size() << ") and number of cells in mesh ("<< (cEnd - cStart) << ").";
            throw std::runtime_error(msg.str());
        } // if
        for (PetscInt c = cStart; c < cEnd; ++c) {
            err = DMSetLabelValue(dmMesh, labelName, c, materialIds[c-cStart]);PYLITH_CHECK_ERROR(err);
        } // for
    } else {
        err = DMCreateLabel(dmMesh, labelName);PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // _setMaterials


// ----------------------------------------------------------------------
// Get material identifiers for cells.
void
pylith::meshio::MeshIO::_getMaterials(int_array* materialIds) const {
    PYLITH_METHOD_BEGIN;

    assert(materialIds);
    assert(_mesh);

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    materialIds->resize(cellsStratum.size());
    PetscErrorCode err = 0;
    PetscInt matId = 0;
    const char* const labelName = pylith::topology::Mesh::cells_label_name;
    for (PetscInt c = cStart, index = 0; c < cEnd; ++c) {
        err = DMGetLabelValue(dmMesh, labelName, c, &matId);PYLITH_CHECK_ERROR(err);
        (*materialIds)[index++] = matId;
    } // for

    PYLITH_METHOD_END;
} // _getMaterials


// ----------------------------------------------------------------------
// Get names of all groups in mesh.
void
pylith::meshio::MeshIO::_getGroupNames(string_vector* names) const { // _getGroups
    PYLITH_METHOD_BEGIN;

    assert(names);
    assert(_mesh);

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    PetscInt numLabels = 0;
    PetscErrorCode err = 0;
    err = DMGetNumLabels(dmMesh, &numLabels);PYLITH_CHECK_ERROR(err);
    const PetscInt numGroups = numLabels - 3; // Remove depth, celltype, and material labels.
    names->resize(numGroups);

    const std::string& materialLabelName = pylith::topology::Mesh::cells_label_name;
    for (int iGroup = 0, iLabel = 0; iLabel < numLabels; ++iLabel) {
        const char* labelName = NULL;
        err = DMGetLabelName(dmMesh, iLabel, &labelName);PYLITH_CHECK_ERROR(err);
        if ((std::string(labelName) != std::string("depth"))
            && (std::string(labelName) != std::string("celltype"))
            && (std::string(labelName) != materialLabelName)) {
            (*names)[iGroup++] = labelName;
        } // if
    } // for

    PYLITH_METHOD_END;
} // _getGroups


// ----------------------------------------------------------------------
// Get group entities
void
pylith::meshio::MeshIO::_getGroup(int_array* points,
                                  pylith::meshio::MeshBuilder::GroupPtType* groupType,
                                  const char *name) const { // _getGroup
    PYLITH_METHOD_BEGIN;

    assert(points);
    assert(groupType);
    assert(_mesh);

    PetscDM dmMesh = _mesh->getDM();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();
    const PetscInt numCells = cellsStratum.size();

    topology::Stratum verticesStratum(dmMesh, topology::Stratum::DEPTH, 0);
    const PetscInt vStart = verticesStratum.begin();
    const PetscInt vEnd = verticesStratum.end();

    PetscIS groupIS = NULL;
    const PetscInt* groupIndices = NULL;
    PetscErrorCode err;
    err = DMGetStratumIS(dmMesh, name, 1, &groupIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);

    PetscInt totalSize;
    err = DMGetStratumSize(dmMesh, name, 1, &totalSize);PYLITH_CHECK_ERROR(err);

    *groupType = pylith::meshio::MeshBuilder::VERTEX;
    if (( totalSize > 0) && (( groupIndices[0] >= cStart) && ( groupIndices[0] < cEnd) )) {
        *groupType = pylith::meshio::MeshBuilder::CELL;
    } // if

    PetscInt offset = 0;
    PetscInt pStart = cStart;
    PetscInt pEnd = cEnd;
    if (pylith::meshio::MeshBuilder::VERTEX == *groupType) {
        offset = numCells;
        pStart = vStart;
        pEnd = vEnd;
    } // if

    // Count number of cells/vertices, filtering out edges and faces
    PetscInt groupSize = 0;
    for (PetscInt p = 0; p < totalSize; ++p) {
        if (( groupIndices[p] >= pStart) && ( groupIndices[p] < pEnd) ) {
            ++groupSize;
        } // if
    } // for

    points->resize(groupSize);
    for (PetscInt p = 0; p < groupSize; ++p) {
        if (( groupIndices[p] >= pStart) && ( groupIndices[p] < pEnd) ) {
            (*points)[p] = groupIndices[p]-offset;
        } // if
    } // for
    err = ISRestoreIndices(groupIS, &groupIndices);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&groupIS);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _getGroup


// End of file
