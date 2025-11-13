// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/MeshOps.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/utils/array.hh" // USES integer_array
#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

#include <algorithm> // USES std::sort, std::find
#include <map> // USES std::map

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class _MeshOps {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static pylith::integer createSubdomainMesh;
                static pylith::integer createLowerDimMesh;
                static pylith::integer createFromPoints;
                static pylith::integer nondimensionalize;
                static pylith::integer checkTopology;
                static pylith::integer checkTopologyGeometry;
                static pylith::integer checkTopologySymmetry;
                static pylith::integer checkTopologySkeleton;
                static pylith::integer checkMaterialLabels;

                static bool isInitialized;
            };
        };
    }
}
pylith::utils::EventLogger pylith::topology::_MeshOps::Events::logger;
pylith::integer pylith::topology::_MeshOps::Events::createSubdomainMesh;
pylith::integer pylith::topology::_MeshOps::Events::createLowerDimMesh;
pylith::integer pylith::topology::_MeshOps::Events::createFromPoints;
pylith::integer pylith::topology::_MeshOps::Events::nondimensionalize;
pylith::integer pylith::topology::_MeshOps::Events::checkTopology;
pylith::integer pylith::topology::_MeshOps::Events::checkTopologyGeometry;
pylith::integer pylith::topology::_MeshOps::Events::checkTopologySymmetry;
pylith::integer pylith::topology::_MeshOps::Events::checkTopologySkeleton;
pylith::integer pylith::topology::_MeshOps::Events::checkMaterialLabels;
bool pylith::topology::_MeshOps::Events::isInitialized = false;

void
pylith::topology::_MeshOps::Events::init(void) {
    if (isInitialized) {
        return;
    } // if

    logger.setClassName("MeshOps");
    logger.initialize();
    createSubdomainMesh = logger.registerEvent("PL:MeshOps:createSubdomainMesh");
    createLowerDimMesh = logger.registerEvent("PL:MeshOps:createLowerDimMesh");
    createFromPoints = logger.registerEvent("PL:MeshOps:createFromPoints");
    nondimensionalize = logger.registerEvent("PL:MeshOps:nondimensionalize");
    checkTopology = logger.registerEvent("PL:MeshOps:checkTopology");
    checkTopologyGeometry = logger.registerEvent("PL:MeshOps:checkTopologyGeometry");
    checkTopologySymmetry = logger.registerEvent("PL:MeshOps:checkTopologySymmetry");
    checkTopologySkeleton = logger.registerEvent("PL:MeshOps:checkTopologySkeleton");
    checkMaterialLabels = logger.registerEvent("PL:MeshOps:checkMaterialLabels");

    isInitialized = true;
}


// ---------------------------------------------------------------------------------------------------------------------
// Create subdomain mesh using label.
std::unique_ptr<pylith::topology::Mesh>
pylith::topology::MeshOps::createSubdomainMesh(const pylith::topology::Mesh& mesh,
                                               const char* labelName,
                                               const int labelValue,
                                               const char* componentName) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createSubdomainMesh);

    assert(labelName);

    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscErrorCode err = PETSC_SUCCESS;

    PetscBool hasLabel = PETSC_FALSE;
    err = DMHasLabel(dmDomain, labelName, &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    /* :TODO: Add creation of pointSF for submesh */
    PetscDMLabel dmLabel = nullptr;
    err = DMGetLabel(dmDomain, labelName, &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    PetscBool hasLabelValue = PETSC_FALSE;
    err = DMLabelHasValue(dmLabel, labelValue, &hasLabelValue);PYLITH_CHECK_ERROR(err);
    int hasLabelValueIntLocal = int(hasLabelValue);
    int hasLabelValueInt = 0;
    err = MPI_Allreduce(&hasLabelValueIntLocal, &hasLabelValueInt, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmDomain));PYLITH_CHECK_ERROR(err);
    if (!hasLabelValueInt) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' with label value '"
            << labelValue << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    PetscDM dmSubdomain = nullptr;
    err = DMPlexFilter(dmDomain, dmLabel, labelValue, PETSC_FALSE, PETSC_FALSE, nullptr, &dmSubdomain);PYLITH_CHECK_ERROR(err);

    pylith::integer maxConeSizeLocal = 0, maxConeSize = 0;
    err = DMPlexGetMaxSizes(dmSubdomain, &maxConeSizeLocal, nullptr);PYLITH_CHECK_ERROR(err);
    err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmSubdomain));PYLITH_CHECK_ERROR(err);

    if (maxConeSize <= 0) {
        err = DMDestroy(&dmSubdomain);PYLITH_CHECK_ERROR(err);
        std::ostringstream msg;
        msg << "Error while creating mesh of subdomain. Subdomain mesh '" << labelName
            << "' with label value " << labelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label name and value.\n";
        throw std::runtime_error(msg.str());
    } // if

    pylith::real lengthScale;
    err = DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmSubdomain, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    std::unique_ptr<pylith::topology::Mesh> submesh(new pylith::topology::Mesh());assert(submesh);
    submesh->setCoordSys(mesh.getCoordSys());
    submesh->setDM(dmSubdomain, componentName);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createSubdomainMesh);
    PYLITH_METHOD_RETURN(submesh);
} // createSubdomainMesh


// ---------------------------------------------------------------------------------------------------------------------
// Create lower dimension mesh using label.
std::unique_ptr<pylith::topology::Mesh>
pylith::topology::MeshOps::createLowerDimMesh(const pylith::topology::Mesh& mesh,
                                              const char* labelName,
                                              const int labelValue,
                                              const char* componentName) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createLowerDimMesh);
    assert(labelName);

    if (mesh.getDimension() < 1) {
        PYLITH_JOURNAL_LOGICERROR("Cannot create submesh for mesh with dimension < 1.");
    } // if

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmDomain = mesh.getDM();assert(dmDomain);
    PetscBool hasLabel = PETSC_FALSE;
    err = DMHasLabel(dmDomain, labelName, &hasLabel);PYLITH_CHECK_ERROR(err);
    if (!hasLabel) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    /* TODO: Add creation of pointSF for submesh */
    PetscDMLabel dmLabel = nullptr;
    err = DMGetLabel(dmDomain, labelName, &dmLabel);PYLITH_CHECK_ERROR(err);assert(dmLabel);
    PetscBool hasLabelValue = PETSC_FALSE;
    err = DMLabelHasValue(dmLabel, labelValue, &hasLabelValue);PYLITH_CHECK_ERROR(err);
    int hasLabelValueIntLocal = int(hasLabelValue);
    int hasLabelValueInt = 0;
    err = MPI_Allreduce(&hasLabelValueIntLocal, &hasLabelValueInt, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmDomain));PYLITH_CHECK_ERROR(err);
    if (!hasLabelValueInt) {
        std::ostringstream msg;
        msg << "Could not find group of points '" << labelName << "' with label value '"
            << labelValue << "' in PETSc DM mesh.";
        throw std::runtime_error(msg.str());
    } // if

    pylith::integer labelHasVertices = 0;
    { // TEMPORARY: Continue to support creating lower dimension meshes using labels with vertices.
        PetscIS labelIS = nullptr;
        const pylith::integer* labelPoints = nullptr;
        pylith::integer numPoints = 0;
        err = DMGetStratumIS(dmDomain, labelName, labelValue, &labelIS);PYLITH_CHECK_ERROR(err);
        pylith::integer labelHasVerticesLocal = 0;
        if (labelIS) {
            err = ISGetIndices(labelIS, &labelPoints);PYLITH_CHECK_ERROR(err);
            err = DMGetStratumSize(dmDomain, labelName, labelValue, &numPoints);PYLITH_CHECK_ERROR(err);

            topology::Stratum verticesStratum(dmDomain, topology::Stratum::DEPTH, 0);
            const pylith::integer vStart = verticesStratum.begin();
            const pylith::integer vEnd = verticesStratum.end();
            for (pylith::integer iPoint = 0; iPoint < numPoints; ++iPoint) {
                if ((labelPoints[iPoint] >= vStart) && (labelPoints[iPoint] < vEnd) ) {
                    labelHasVerticesLocal = 1;
                    break;
                } // if
            } // if
            err = ISRestoreIndices(labelIS, &labelPoints);PYLITH_CHECK_ERROR(err);
        } // if
        err = ISDestroy(&labelIS);PYLITH_CHECK_ERROR(err);
        err = MPI_Allreduce(&labelHasVerticesLocal, &labelHasVertices, 1, MPI_INT, MPI_MAX,
                            PetscObjectComm((PetscObject) dmDomain));PYLITH_CHECK_ERROR(err);

        if (labelHasVertices) {
            pythia::journal::warning_t warning("deprecated");
            warning << pythia::journal::at(__HERE__)
                    << "DEPRECATION: Creating lower dimension mesh from label with vertices. "
                    << "This feature will be removed in v6.0. "
                    << "In the future, you will need to mark boundaries not vertices for boundary conditions."
                    << pythia::journal::endl;
        } // if
    } // TEMPORARY

    // We use DMPlexCreateSubmesh() instead of DMPlexFilter, because we want the submesh to have
    // domain cells hanging off of it, which allows us to project from the submesh to the domain mesh
    // to set boundary conditions using the auxiliary fields defined over the submesh.
    // DMPlexCreateSubmesh() requires a completed label.
    PetscDMLabel dmLabelFull = nullptr;
    err = DMLabelDuplicate(dmLabel, &dmLabelFull);PYLITH_CHECK_ERROR(err);
    err = DMPlexLabelComplete(dmDomain, dmLabelFull);PYLITH_CHECK_ERROR(err);

    PetscDM dmSubmesh = nullptr;
    const PetscBool markedFaces = !labelHasVertices ? PETSC_TRUE : PETSC_FALSE;
    err = DMPlexCreateSubmesh(dmDomain, dmLabelFull, labelValue, markedFaces, &dmSubmesh);PYLITH_CHECK_ERROR(err);
    err = DMLabelDestroy(&dmLabelFull);PYLITH_CHECK_ERROR(err);

    pylith::integer maxConeSizeLocal = 0, maxConeSize = 0;
    err = DMPlexGetMaxSizes(dmSubmesh, &maxConeSizeLocal, nullptr);PYLITH_CHECK_ERROR(err);
    err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                        PetscObjectComm((PetscObject) dmSubmesh));PYLITH_CHECK_ERROR(err);

    if (maxConeSize <= 0) {
        err = DMDestroy(&dmSubmesh);PYLITH_CHECK_ERROR(err);
        std::ostringstream msg;
        msg << "Error while creating lower dimension mesh. Submesh '" << labelName
            << "' with label value " << labelValue << " does not contain any cells.\n"
            << "Check that you are using the correct label name and value.\n";
        throw std::runtime_error(msg.str());
    } // if

    // Set length scale
    pylith::real lengthScale;
    err = DMPlexGetScale(dmDomain, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmSubmesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    std::unique_ptr<pylith::topology::Mesh> submesh(new pylith::topology::Mesh());assert(submesh);
    submesh->setCoordSys(mesh.getCoordSys());
    submesh->setDM(dmSubmesh, componentName);

    // Check topology
    MeshOps::checkTopology(*submesh);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createLowerDimMesh);
    PYLITH_METHOD_RETURN(submesh);
} // createLowerDimMesh


// ---------------------------------------------------------------------------------------------------------------------
// Create 0-dimension mesh from points.
std::unique_ptr<pylith::topology::Mesh>
pylith::topology::MeshOps::createFromPoints(const pylith::real* points,
                                            const size_t numPoints,
                                            const std::shared_ptr<spatialdata::geocoords::CoordSys>& cs,
                                            const pylith::real lengthScale,
                                            MPI_Comm comm,
                                            const char* componentName) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::createFromPoints);
    assert(cs);

    PetscErrorCode err = PETSC_SUCCESS;

    const int meshDim = 0;
    std::unique_ptr<pylith::topology::Mesh> mesh(new pylith::topology::Mesh(meshDim, comm));assert(mesh);

    PetscDM dmPoints = nullptr;
    const pylith::integer depth = 0;
    pylith::integer dmNumPoints[1];
    dmNumPoints[0] = numPoints;
    pylith::integer_array dmConeSizes(0, numPoints);
    pylith::integer_array dmCones(0, numPoints);
    pylith::integer_array dmConeOrientations(0, numPoints);

    const size_t spaceDim = cs->getSpaceDim();

    err = DMPlexCreate(comm, &dmPoints);PYLITH_CHECK_ERROR(err);
    err = DMSetDimension(dmPoints, 0);PYLITH_CHECK_ERROR(err);
    err = DMSetCoordinateDim(dmPoints, spaceDim);PYLITH_CHECK_ERROR(err);
    if (numPoints > 0) {
        err = DMPlexCreateFromDAG(dmPoints, depth, dmNumPoints, &dmConeSizes[0], &dmCones[0],
                                  &dmConeOrientations[0], points);PYLITH_CHECK_ERROR(err);
    } else {
        pylith::integer empty[1];
        empty[0] = 0;
        err = DMPlexCreateFromDAG(dmPoints, depth, dmNumPoints, &empty[0], &empty[0],
                                  &empty[0], points);PYLITH_CHECK_ERROR(err);
    } // if/else

    PetscSF sf = nullptr;
    err = DMGetPointSF(dmPoints, &sf);PYLITH_CHECK_ERROR(err);
    err = PetscSFSetGraph(sf, numPoints, 0, nullptr, PETSC_COPY_VALUES, nullptr, PETSC_COPY_VALUES);

    mesh->setCoordSys(cs);
    mesh->setDM(dmPoints, componentName);
    err = DMPlexSetScale(mesh->getDM(), PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::createFromPoints);
    PYLITH_METHOD_RETURN(mesh);
} // createFromPoints


// ---------------------------------------------------------------------------------------------------------------------
// Nondimensionalize the finite-element mesh.
void
pylith::topology::MeshOps::nondimensionalize(Mesh* const mesh,
                                             const spatialdata::units::Nondimensional& normalizer) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::nondimensionalize);

    assert(mesh);

    PetscVec coordVec = nullptr;
    const pylith::real lengthScale = normalizer.getLengthScale();
    PetscErrorCode err = PETSC_SUCCESS;

    PetscDM dmMesh = mesh->getDM();assert(dmMesh);
    err = DMGetCoordinatesLocal(dmMesh, &coordVec);PYLITH_CHECK_ERROR(err);assert(coordVec);
    err = VecScale(coordVec, 1.0/lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmMesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dmMesh, nullptr, "-pylith_nondim_dm_view");PYLITH_CHECK_ERROR(err);

    const pylith::integer dim = mesh->getDimension();
    if (dim < 1) {
        PYLITH_METHOD_END;
    } // if
    pylith::real coordMin[3];
    pylith::real coordMax[3];
    err = DMGetBoundingBox(dmMesh, coordMin, coordMax);
    pylith::real volume = 1.0;
    for (int i = 0; i < dim; ++i) {
        volume *= coordMax[i] - coordMin[i];
    } // for
    assert(dim > 0);
    const pylith::real avgCellDim = pow(volume / MeshOps::getNumCells(*mesh), 1.0/dim);
    const pylith::real avgDimTolerance = 0.02;
    if (avgCellDim < avgDimTolerance) {
        std::ostringstream msg;
        msg << "Nondimensional average cell dimension (" << avgCellDim << ") is less than minimum tolerance ("
            << avgDimTolerance << "). This usually means the length scale (" << lengthScale << ") used in the "
            << "nondimensionalization needs to be smaller. Based on the average cell size, a value of about "
            << pow(10, int(log10(avgCellDim*lengthScale))) << " should be appropriate.";
        throw std::runtime_error(msg.str());
    } // if/else

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::nondimensionalize);
    PYLITH_METHOD_END;
} // nondimensionalize


// ---------------------------------------------------------------------------------------------------------------------
// Strip out "ghost" cells hanging off mesh
PetscDM
pylith::topology::MeshOps::removeHangingCells(const PetscDM& dmMesh) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = PETSC_SUCCESS;
    PetscDM dmClean = PETSC_nullptrPTR;

    MPI_Comm comm = PetscObjectComm((PetscObject) dmMesh);
    pylith::topology::Stratum cells(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    DMPolytopeType cellType;
    err = DMPlexGetCellType(dmMesh, cells.begin(), &cellType);PYLITH_CHECK_ERROR(err);
    if (DMPolytopeTypeGetDim(cellType) < 0) {
        // Hanging cells have dim == -1

        // Create label over cells 1 dimension lower
        PetscDMLabel labelInclude = PETSC_nullptrPTR;
        const pylith::integer labelValue = 1;
        err = DMLabelCreate(comm, "no_hanging_cells", &labelInclude);PYLITH_CHECK_ERROR(err);
        pylith::topology::Stratum faces(dmMesh, pylith::topology::Stratum::HEIGHT, 1);
        for (pylith::integer face = faces.begin(); face < faces.end(); ++face) {
            err = DMLabelSetValue(labelInclude, face, labelValue);PYLITH_CHECK_ERROR(err);
        } // for

        err = DMPlexFilter(dmMesh, labelInclude, labelValue, PETSC_FALSE, PETSC_FALSE, PETSC_nullptrPTR, &dmClean);PYLITH_CHECK_ERROR(err);
        err = DMLabelDestroy(&labelInclude);PYLITH_CHECK_ERROR(err);

        // Create section using subpoint map to ensure sections are consistent.
        PetscIS subpointIS = PETSC_nullptrPTR;
        PetscSection sectionOld = PETSC_nullptrPTR, sectionNew = PETSC_nullptrPTR;
        err = DMPlexGetSubpointIS(dmClean, &subpointIS);PYLITH_CHECK_ERROR(err);
        err = DMGetLocalSection(dmMesh, &sectionOld);PYLITH_CHECK_ERROR(err);
        err = PetscSectionCreateSubmeshSection(sectionOld, subpointIS, &sectionNew);PYLITH_CHECK_ERROR(err);
        err = DMSetLocalSection(dmClean, sectionNew);PYLITH_CHECK_ERROR(err);
        err = PetscSectionDestroy(&sectionNew);PYLITH_CHECK_ERROR(err);

        pylith::real lengthScale = 0.0;
        err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
        err = DMPlexSetScale(dmClean, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
    } else {
        dmClean = dmMesh;
        err = PetscObjectReference((PetscObject) dmClean);
    } // if/else

    PYLITH_METHOD_RETURN(dmClean);
}


// ---------------------------------------------------------------------------------------------------------------------
// Check topology of mesh.
void
pylith::topology::MeshOps::checkTopology(const Mesh& mesh) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopology);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);

    DMLabel subpointMap;
    PetscErrorCode ierr = DMPlexGetSubpointMap(dmMesh, &subpointMap);PYLITH_CHECK_ERROR(ierr);
    pylith::integer cellHeight = subpointMap ? 1 : 0;

    PetscErrorCode err = PETSC_SUCCESS;
    err = DMViewFromOptions(dmMesh, nullptr, "-pylith_checktopo_dm_view");PYLITH_CHECK_ERROR(err);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologyGeometry);
    err = DMPlexCheckGeometry(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of the mesh.");
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologyGeometry);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologySymmetry);
    err = DMPlexCheckSymmetry(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of mesh associated with symmetry of adjacency information.");
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologySymmetry);

    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkTopologySkeleton);
    err = DMPlexCheckSkeleton(dmMesh, cellHeight);PYLITH_CHECK_ERROR_MSG(err, "Error in topology of mesh cells.");
    err = DMPlexCheckOrphanVertices(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Mesh contains vertices not connected to cells.");
    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopologySkeleton);

    /* Other check functions that we are not using:
     *
     * DMPlexCheckFaces() - not compatible with cohesive cells.
     *
     * DMPlexCheckInterfaceCones() - very slow
     */

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkTopology);
    PYLITH_METHOD_END;
} // checkTopology


// ---------------------------------------------------------------------------------------------------------------------
bool
pylith::topology::MeshOps::isSimplexMesh(const Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    bool isSimplex = false;

    PetscErrorCode err = PETSC_SUCCESS;
    const PetscDM dm = mesh.getDM();
    pylith::integer vStart = 0, vEnd = 0;
    err = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    if (vStart != vEnd) { // Test for simplex only works if we have points.
        pylith::integer closureSize = 0;
        pylith::integer* closure = nullptr;
        const int dim = mesh.getDimension();

        err = DMPlexGetTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        pylith::integer numVertices = 0;
        for (pylith::integer c = 0; c < closureSize*2; c += 2) {
            if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
                ++numVertices;
            } // if
        } // for
        if (numVertices == dim+1) {
            isSimplex = PETSC_TRUE;
        } // if
        err = DMPlexRestoreTransitiveClosure(dm, 0, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // if

    // Communicate result of isSimplex to all processes.
    int intSimplexLocal = isSimplex ? 1 : 0;
    int intSimplexGlobal = 0;
    MPI_Allreduce(&intSimplexLocal, &intSimplexGlobal, 1, MPI_INT, MPI_LOR, mesh.getComm());
    isSimplex = intSimplexGlobal == 1;

    PYLITH_METHOD_RETURN(isSimplex);
} // isSimplexMesh


// ---------------------------------------------------------------------------------------------------------------------
bool
pylith::topology::MeshOps::isCohesiveCell(const PetscDM dm,
                                          const pylith::integer cell) {
    bool isCohesive = false;

    DMPolytopeType ct;
    PetscErrorCode err = DMPlexGetCellType(dm, cell, &ct);PYLITH_CHECK_ERROR(err);
    if ((ct == DM_POLYTOPE_SEG_PRISM_TENSOR) ||
        (ct == DM_POLYTOPE_TRI_PRISM_TENSOR) ||
        (ct == DM_POLYTOPE_QUAD_PRISM_TENSOR)) { isCohesive = true; }

    return isCohesive;
} // isCohesiveCell


// ----------------------------------------------------------------------
// Get number of vertices in mesh.
pylith::integer
pylith::topology::MeshOps::getNumVertices(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::integer nvertices = 0;
    pylith::integer begin = 0, end = 0;
    PetscErrorCode err = DMPlexGetDepthStratum(dmMesh, 0, &begin, &end);PYLITH_CHECK_ERROR(err);
    nvertices = end-begin;

    PYLITH_METHOD_RETURN(nvertices);
}


// ----------------------------------------------------------------------
// Get number of cells in mesh.
pylith::integer
pylith::topology::MeshOps::getNumCells(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    pylith::integer ncells = 0;
    pylith::integer begin = 0, end = 0;
    const int cellHeight = 0;
    PetscErrorCode err = DMPlexGetHeightStratum(dmMesh, cellHeight, &begin, &end);PYLITH_CHECK_ERROR(err);
    ncells = end-begin;

    PYLITH_METHOD_RETURN(ncells);
}


// ----------------------------------------------------------------------
// Get number of vertices in a cell.
pylith::integer
pylith::topology::MeshOps::getNumCorners(const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    pylith::integer numCorners = 0;
    PetscDM dmMesh = mesh.getDM();assert(dmMesh);

    pylith::integer cStart, cEnd, vStart, vEnd, closureSize, *closure = nullptr;
    PetscErrorCode err = PETSC_SUCCESS;
    const int cellHeight = 0;
    err = DMPlexGetHeightStratum(dmMesh, cellHeight, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    if (cEnd > cStart) {
        err = DMPlexGetTransitiveClosure(dmMesh, cStart, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
        for (pylith::integer c = 0; c < closureSize*2; c += 2) {
            if ((closure[c] >= vStart) && (closure[c] < vEnd)) {++numCorners;}
        } // for
        err = DMPlexRestoreTransitiveClosure(dmMesh, cStart, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(numCorners);
}


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::topology::MeshOps::checkMaterialLabels(const pylith::topology::Mesh& mesh,
                                               pylith::integer_array& labelValues) {
    PYLITH_METHOD_BEGIN;
    _MeshOps::Events::init();
    _MeshOps::Events::logger.eventBegin(_MeshOps::Events::checkMaterialLabels);

    PetscErrorCode err = PETSC_SUCCESS;

    // Create map with indices for each material
    const size_t numIds = labelValues.size();
    std::map<int, int> materialIndex;
    for (size_t i = 0; i < numIds; ++i) {
        materialIndex[labelValues[i]] = i;
    } // for

    integer_array matCellCounts(numIds);
    matCellCounts = 0;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    Stratum cellsStratum(dmMesh, Stratum::HEIGHT, 0);
    const pylith::integer cStart = cellsStratum.begin();
    const pylith::integer cEnd = cellsStratum.end();

    PetscDMLabel materialsLabel = nullptr;
    const char* const labelName = pylith::topology::Mesh::cells_label_name;
    err = DMGetLabel(dmMesh, labelName, &materialsLabel);PYLITH_CHECK_ERROR(err);assert(materialsLabel);

    int *matBegin = &labelValues[0];
    int *matEnd = &labelValues[0] + labelValues.size();
    std::sort(matBegin, matEnd);

    for (pylith::integer c = cStart; c < cEnd; ++c) {
        pylith::integer matId;

        err = DMLabelGetValue(materialsLabel, c, &matId);PYLITH_CHECK_ERROR(err);
        if (matId < 0) {
            // :KLUDGE: Skip cells that are probably hybrid cells in halo
            // around fault that we currently ignore when looping over
            // materials (including cohesive cells).
            continue;
        } // if
        const int *result = std::find(matBegin, matEnd, matId);
        if (result == matEnd) {
            std::ostringstream msg;
            msg << "Material label_value '" << matId << "' for cell '" << c
                << "' does not match the label_value of any materials or interfaces.";
            throw std::runtime_error(msg.str());
        } // if

        const size_t matIndex = materialIndex[matId];
        assert(0 <= matIndex && matIndex < numIds);
        ++matCellCounts[matIndex];
    } // for

    // Make sure each material has cells.
    integer_array matCellCountsAll(matCellCounts.size());
    err = MPI_Allreduce(&matCellCounts[0], &matCellCountsAll[0],
                        matCellCounts.size(), MPI_INT, MPI_SUM, mesh.getComm());PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < numIds; ++i) {
        const int matId = labelValues[i];
        const size_t matIndex = materialIndex[matId];
        assert(0 <= matIndex && matIndex < numIds);
        if (matCellCountsAll[matIndex] <= 0) {
            std::ostringstream msg;
            msg << "No cells associated with material with id '" << matId << "'.";
            throw std::runtime_error(msg.str());
        } // if
    } // for

    _MeshOps::Events::logger.eventEnd(_MeshOps::Events::checkMaterialLabels);
    PYLITH_METHOD_END;
} // checkMaterialIds


// End of file
