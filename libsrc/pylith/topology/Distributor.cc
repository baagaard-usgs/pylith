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

#include "pylith/topology/Distributor.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/meshio/OutputSubfield.hh" // USES OutputSubfield
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/faults/FaultCohesive.hh" // USES FaultCohesive
#include "pylith/meshio/DataWriter.hh" // USES DataWriter
#include "pylith/utils/journals.hh" // pythia::journal

#include <cstring> // USES strlen()
#include <strings.h> // USES strcasecmp()
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class _Distributor {
public:

            static
            void write(pylith::meshio::DataWriter* writer,
                       const pylith::topology::Mesh& mesh);

            /** Distribute custom overlap based on PETSc labels.
             *
             * The overlap excludes cohesive cells but includes cells adjacent to faults.
             * This is a custom version of DMPlexDistributeOverlap()
             *
             * @param[out] dmOverlap PETSc DM for the overlap.
             * @param[in] dmMesh PETSc DM for the current mesh.
             * @param[in] faults Array of fault interfaces.
             *
             * @returns PETSc error code (0==success).
             */
            static
            PetscErrorCode distributeOverlap(PetscDM* dmOverlap,
                                             PetscDM dmMesh,
                                             const std::vector<pylith::faults::FaultCohesive*>& faults);

        }; // _Distributor
    } // topology
} // pylith

// ------------------------------------------------------------------------------------------------
static PetscErrorCode
DMPlexGetAdjacency_SupportOnly_Internal(DM dm,
                                        PetscInt p,
                                        PetscInt *adjSize,
                                        PetscInt adj[],
                                        void *ctx) {
    const PetscInt *support = NULL;
    PetscInt maxAdjSize = *adjSize;

    PetscFunctionBeginHot;
    PylithCallPetsc(DMPlexGetSupportSize(dm, p, adjSize));
    PylithCallPetsc(DMPlexGetSupport(dm, p, &support));
    PetscCheck(*adjSize <= maxAdjSize, PETSC_COMM_SELF, PETSC_ERR_PLIB, "Invalid mesh exceeded adjacency allocation (%" PetscInt_FMT ")", maxAdjSize);
    for (PetscInt s = 0; s < *adjSize; ++s) {
        adj[s] = support[s];
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::Distributor::Distributor(void) :
    _writer(nullptr),
    _partitioner("parmetis"),
    _useEdgeWeighting(true) {}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::Distributor::~Distributor(void) {
    _writer = nullptr; // :KLUDGE: Use shared pointer
}


// ------------------------------------------------------------------------------------------------
// Set data writer.
void
pylith::topology::Distributor::setDataWriter(pylith::meshio::DataWriter* writer) {
    _writer = writer;
} // setDataWriter


// ------------------------------------------------------------------------------------------------
// Set edge weighting.
void
pylith::topology::Distributor::setUseEdgeWeighting(const bool flag) {
    _useEdgeWeighting = flag;
} // setUseEdgeWeighting


// ------------------------------------------------------------------------------------------------
// Set partitioner.
void
pylith::topology::Distributor::setPartitioner(const char* partitioner) {
    if ((0 == strcasecmp(partitioner, "parmetis")) || (0 == strcasecmp(partitioner, "chaco")) || (0 == strcasecmp(partitioner, "simple"))) {
        _partitioner = partitioner;
    } else {
        std::ostringstream msg;
        msg << "Unknown partitioner '" << partitioner << "'. Partitioner must be 'parmetis', 'chaco', or 'simple'.";
        throw std::runtime_error(msg.str());
    } // if/else
} // setPartitioner


// ------------------------------------------------------------------------------------------------
// Distribute mesh.
pylith::topology::Mesh*
pylith::topology::Distributor::distribute(const pylith::topology::Mesh& mesh,
                                          const std::vector<pylith::faults::FaultCohesive*>& faults) const {
    PYLITH_METHOD_BEGIN;

    pylith::topology::Mesh* newMesh = new pylith::topology::Mesh();assert(newMesh);
    newMesh->setCoordSys(mesh.getCoordSys());

    const int commRank = mesh.getCommRank();
    if (0 == commRank) {
        PYLITH_COMPONENT_INFO("Partitioning mesh using PETSc '" << _partitioner << "' partitioner.");
    } // if

    PetscPartitioner partitioner = nullptr;
    PetscDM dmOrig = mesh.getDM();assert(dmOrig);
    PylithCallPetsc(DMPlexGetPartitioner(dmOrig, &partitioner));
    PylithCallPetsc(PetscPartitionerSetType(partitioner, _partitioner.c_str()));
    if ((_partitioner == std::string("parmetis")) && _useEdgeWeighting) {
        PylithCallPetsc(PetscOptionsSetValue(NULL, "-petscpartitioner_use_vertex_weights", "true"));
        PylithCallPetsc(PetscPartitionerSetFromOptions(partitioner));
    } // if

    if (0 == commRank) {
        PYLITH_COMPONENT_INFO("Distributing partitioner mesh.");
    } // if

    PetscDM dmTmp = NULL, dmNew = NULL;
    const PetscInt overlap = 0;
    PylithCallPetsc(DMPlexDistribute(dmOrig, overlap, NULL, &dmTmp));
    PylithCallPetsc(_Distributor::distributeOverlap(&dmNew, dmTmp, faults));
    PylithCallPetsc(DMDestroy(&dmTmp));
    if (dmNew) {
        PylithCallPetsc(DMPlexDistributeSetDefault(dmNew, PETSC_FALSE));
        PylithCallPetsc(DMPlexReorderCohesiveSupports(dmNew));
        PylithCallPetsc(DMViewFromOptions(dmNew, NULL, "-pylith_dist_dm_view"));
        newMesh->setDM(dmNew, "domain");
    } else {
        PetscObjectReference(PetscObject(dmOrig));
        newMesh->setDM(dmOrig);
    } // if/else

    pythia::journal::debug_t debug(PyreComponent::getName());
    if (debug.state()) {
        newMesh->view(":mesh_distributed.txt:ascii_info_detail");
        newMesh->view(":mesh_distributed.tex:ascii_latex");
    } // if
    if (_writer) {
        _Distributor::write(_writer, *newMesh);
    } // if

    PYLITH_METHOD_RETURN(newMesh);
} // distribute


// ------------------------------------------------------------------------------------------------
// Write partitioning info for distributed mesh.
void
pylith::topology::_Distributor::write(meshio::DataWriter* const writer,
                                      const topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    // Setup and allocate PETSc vector
    const int commRank = mesh.getCommRank();
    PylithScalar rankReal = PylithReal(commRank);

    pylith::topology::Field partitionField(mesh);
    partitionField.setLabel("partition");

    pylith::topology::Field::Description description;
    description.label = "partition";
    description.alias = "partition";
    description.vectorFieldType = pylith::topology::Field::SCALAR;
    description.numComponents = 1;
    description.componentNames.resize(1);
    description.componentNames[0] = "rank";
    description.scale = 1.0;
    description.validator = NULL;

    pylith::topology::Field::Discretization discretization(0, 1);

    partitionField.subfieldAdd(description, discretization);
    partitionField.subfieldsSetup();
    partitionField.createDiscretization();
    partitionField.allocate();
    partitionField.zeroLocal();
    partitionField.createOutputVector();

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    topology::Stratum cellsStratum(dmMesh, pylith::topology::Stratum::HEIGHT, 0);
    const PetscInt cStart = cellsStratum.begin();
    const PetscInt cEnd = cellsStratum.end();

    VecVisitorMesh partitionVisitor(partitionField);
    PetscScalar* partitionArray = partitionVisitor.localArray();
    for (PetscInt point = cStart; point < cEnd; ++point) {
        const PetscInt off = partitionVisitor.sectionOffset(point);
        if (partitionVisitor.sectionDof(point) > 0) {
            partitionArray[off] = rankReal;
        } // if
    } // for
    partitionVisitor.clear();
    partitionField.scatterLocalToOutput();

    const int basisOrder = 0;
    const int refineLevels = 0;
    pylith::meshio::OutputSubfield* outputField =
        pylith::meshio::OutputSubfield::create(partitionField, mesh, "partition", basisOrder, refineLevels);
    outputField->project(partitionField.getOutputVector());

    const PylithScalar t = 0.0;
    const bool isInfo = true;
    writer->open(mesh, isInfo);
    writer->openTimeStep(t, mesh);
    writer->writeCellField(t, *outputField);
    writer->closeTimeStep();
    writer->close();

    delete outputField;outputField = NULL;

    PYLITH_METHOD_END;
} // write


// ------------------------------------------------------------------------------------------------
// This is a copy of DMPlexDistributeOverlap()
PetscErrorCode
pylith::topology::_Distributor::distributeOverlap(PetscDM* dmOverlap,
                                                  PetscDM dmMesh,
                                                  const std::vector<pylith::faults::FaultCohesive*>& faults) {
    PYLITH_METHOD_BEGIN;
    assert(dmOverlap);

    MPI_Comm comm;
    PetscMPIInt size, rank;
    PetscSection rootSection, leafSection;
    PetscIS rootrank, leafrank;
    PetscDM dmCoord;
    PetscDMLabel lblOverlap;
    PetscSF sfOverlap, sfStratified, sfPoint;
    PetscInt dim;

    const size_t numFaults = faults.size();
    if (0 == numFaults) {
        PylithCallPetsc(PetscObjectReference((PetscObject)dmMesh));
        *dmOverlap = dmMesh;
        PYLITH_METHOD_RETURN(0);
    } // if

    PetscDMLabel* ovIncludeLabels = (numFaults > 0) ? new PetscDMLabel[numFaults] : NULL;
    PetscInt* ovIncludeLabelValues = (numFaults > 0) ? new PetscInt[numFaults] : NULL;
    PetscDMLabel* ovExcludeLabels = (numFaults > 0) ? new PetscDMLabel[numFaults] : NULL;
    PetscInt* ovExcludeLabelValues = (numFaults > 0) ? new PetscInt[numFaults] : NULL;

    PylithCallPetsc(DMGetDimension(dmMesh, &dim));
    for (int i = 0; i < numFaults; ++i) {
        const char* surfaceLabelName = faults[i]->getSurfaceLabelName();
        PylithCallPetsc(DMGetLabel(dmMesh, surfaceLabelName, &ovIncludeLabels[i]));
        ovIncludeLabelValues[i] = dim - 1; // faults[i]->getSurfaceLabelValue();

        const char* cohesiveLabelName = faults[i]->getCohesiveLabelName();
        PylithCallPetsc(DMGetLabel(dmMesh, cohesiveLabelName, &ovExcludeLabels[i]));
        ovExcludeLabelValues[i] = faults[i]->getCohesiveLabelValue();
    } // for

    PylithCallPetsc(PetscObjectGetComm((PetscObject)dmMesh,&comm));
    PetscCallMPI(MPI_Comm_size(comm, &size));
    PetscCallMPI(MPI_Comm_rank(comm, &rank));
    /* Compute point overlap with neighbouring processes on the distributed DM */
    PylithCallPetsc(PetscSectionCreate(comm, &rootSection));
    PylithCallPetsc(PetscSectionCreate(comm, &leafSection));
    PylithCallPetsc(DMPlexSetAdjacencyUser(dmMesh, DMPlexGetAdjacency_SupportOnly_Internal, NULL));
    PylithCallPetsc(DMPlexDistributeOwnership(dmMesh, rootSection, &rootrank, leafSection, &leafrank));
    PylithCallPetsc(DMPlexCreateOverlapLabelFromLabels(dmMesh, numFaults, ovIncludeLabels, ovIncludeLabelValues,
                                                       numFaults, ovExcludeLabels, ovExcludeLabelValues, rootSection, rootrank, leafSection, leafrank, &lblOverlap));
    PylithCallPetsc(DMPlexSetAdjacencyUser(dmMesh, NULL, NULL));

    delete[] ovIncludeLabels;ovIncludeLabels = NULL;
    delete[] ovIncludeLabelValues;ovIncludeLabelValues = NULL;
    delete[] ovExcludeLabels;ovExcludeLabels = NULL;
    delete[] ovExcludeLabelValues;ovExcludeLabelValues = NULL;

    /* Convert overlap label to stratified migration SF */
    PylithCallPetsc(DMPlexPartitionLabelCreateSF(dmMesh, lblOverlap, PETSC_TRUE, &sfOverlap));
    PylithCallPetsc(DMPlexStratifyMigrationSF(dmMesh, sfOverlap, &sfStratified));
    PylithCallPetsc(PetscSFDestroy(&sfOverlap));
    sfOverlap = sfStratified;
    PylithCallPetsc(PetscObjectSetName((PetscObject) sfOverlap, "Overlap SF"));
    PylithCallPetsc(PetscSFSetFromOptions(sfOverlap));

    PylithCallPetsc(PetscSectionDestroy(&rootSection));
    PylithCallPetsc(PetscSectionDestroy(&leafSection));
    PylithCallPetsc(ISDestroy(&rootrank));
    PylithCallPetsc(ISDestroy(&leafrank));

    /* Build the overlapping DM */
    PylithCallPetsc(DMPlexCreate(comm, dmOverlap));
    PylithCallPetsc(PetscObjectSetName((PetscObject) *dmOverlap, "Parallel Mesh"));
    PylithCallPetsc(DMPlexMigrate(dmMesh, sfOverlap, *dmOverlap));
    /* Store the overlap in the new DM */
    PylithCallPetsc(DMPlexSetOverlap(*dmOverlap, dmMesh, 1));
    /* Build the new point SF */
    PylithCallPetsc(DMPlexCreatePointSF(*dmOverlap, sfOverlap, PETSC_FALSE, &sfPoint));
    PylithCallPetsc(DMSetPointSF(*dmOverlap, sfPoint));
    PylithCallPetsc(DMGetCoordinateDM(*dmOverlap, &dmCoord));
    if (dmCoord) { PylithCallPetsc(DMSetPointSF(dmCoord, sfPoint));}
    PylithCallPetsc(PetscSFDestroy(&sfPoint));
    /* Cleanup overlap partition */
    PylithCallPetsc(DMLabelDestroy(&lblOverlap));
    PylithCallPetsc(PetscSFDestroy(&sfOverlap));

    PYLITH_METHOD_RETURN(0);
} // distributeOverlap


// End of file
