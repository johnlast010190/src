/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.
    And the precice adapter solver object is based on the preCICE-
    adapter for OpenFOAM.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2017-2023 Gerasimos Chourdakis
    (c) 2023 Esi Ltd.

Class
    preciceAdapter::Adapter

SourceFiles
    Adapter.C

\*---------------------------------------------------------------------------*/

#ifndef PRECICEADAPTER_H
#define PRECICEADAPTER_H

#include "Interface.H"

// Conjugate Heat Transfer module
#include "CHT/CHT.H"

// Fluid-Structure Interaction module
#include "FSI/FSI.H"

// Fluid-Fluid module
#include "FF/FF.H"

// NOTE: If you want to couple a new variable, include your module's header here.
// You also need to include it in the Make/files file.
// In case you use additional OpenFOAM symbols, you may also need to specify
// the respective libraries in the Make/options.

// OpenFOAM header files
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"

// preCICE Solver Interface
#include "SolverInterface.hpp"

namespace preciceAdapter
{

class Adapter
{

private:
    //- Structure of the configuration of each coupling interface.
    //  Every interface needs to know the coupling mesh, the OpenFOAM
    //  patches that consist the coupling surface and the kinds
    //  of data that are exchanged.
    struct InterfaceConfig
    {
        std::string meshName;
        std::string locationsType;
        bool meshConnectivity;
        std::vector<std::string> patchNames;
        std::vector<std::string> writeData;
        std::vector<std::string> readData;
    };

    //- Configuration interfaces
    std::vector<struct InterfaceConfig> interfacesConfig_;

#ifdef ADAPTER_ENABLE_TIMINGS
    //- Timers
    clockValue timeInConfigRead_;
    clockValue timeInPreciceConstruct_;
    clockValue timeInMeshSetup_;
    clockValue timeInInitialize_;
    clockValue timeInInitializeData_;
    clockValue timeInCheckpointingSetup_;
    clockValue timeInWrite_;
    clockValue timeInAdvance_;
    clockValue timeInRead_;
    clockValue timeInCheckpointingRead_;
    clockValue timeInCheckpointingWrite_;
    clockValue timeInWriteResults_;
    clockValue timeInFinalize_;
#endif

    //- OpenFOAM runTime object
    const Foam::Time& runTime_;

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    // Configuration parameters used in the Adapter

    //- preCICE configuration file name
    Foam::fileName preciceConfigFilename_;

    //- preCICE participant name
    Foam::word participantName_;

    //- Switch for checkpointing
    bool checkpointing_ = false;

    //- Switch to enable the ConjugateHeatTransfer module
    bool CHTenabled_ = false;

    //- Switch to enable the FluidStructureInteraction module
    bool FSIenabled_ = false;

    //- Switch to enable the FluidFluid module
    bool FFenabled_ = false;

    // NOTE: Add a switch for your new module here

    //- Interfaces
    std::vector<Interface*> interfaces_;

    //- preCICE solver interface
    precice::SolverInterface* precice_ = NULL;

    //- preCICE solver interface initialized
    bool preciceInitialized_ = false;

    //- Conjugate Heat Transfer module object
    CHT::ConjugateHeatTransfer* CHT_ = NULL;

    //- Fluid-Structure Interaction module object
    FSI::FluidStructureInteraction* FSI_ = NULL;

    //- Fluid-Fluid module object
    FF::FluidFluid* FF_ = NULL;

    // NOTE: Add here a pointer for your new module object

    // Timesteps

    //- Timestep dictated by preCICE
    double timestepPrecice_;

    // Checkpointing

    //- Checkpointed time (value)
    Foam::scalar couplingIterationTimeValue_;

    //- Checkpointed time (index)
    Foam::label couplingIterationTimeIndex_;

    //- Checkpointed mesh points
    Foam::pointField meshPoints_;
    Foam::pointField oldMeshPoints_;
    bool meshCheckPointed = false;

    // TODO: Currently unused, see storeMeshPoints().
    //- Checkpointed mesh volume
    // bool oldVolsStored = false;
    // Foam::volScalarField::Internal * oldVols_;
    // Foam::volScalarField::Internal * oldOldVols_;
    // int curTimeIndex_ = 0;

    // Vectors of pointers to the checkpointed mesh fields and their copies

    //- Checkpointed surfaceScalarField mesh fields (copies)
    Foam::PtrList<Foam::surfaceScalarField> meshSurfaceScalarFieldCopies_;

    //- Checkpointed surfaceVectorField mesh fields (copies)
    Foam::PtrList<Foam::surfaceVectorField> meshSurfaceVectorFieldCopies_;

    //- Checkpointed volVectorField mesh fields (copies)
    Foam::PtrList<Foam::volVectorField> meshVolVectorFieldCopies_;

    // TODO checkpoint for the V0 (Old volume) and V00 (Old-Old volume) fields.
    //- Checkpointed volScalarField mesh fields
    Foam::PtrList<Foam::volScalarField::Internal> volScalarInternalFields_;

    //- Checkpointed volScalarField mesh fields (copies)
    Foam::PtrList<Foam::volScalarField::Internal> volScalarInternalFieldCopies_;

    // Vectors of pointers to the checkpointed fields and their copies

    //- Checkpointed volScalarField fields (copies)
    Foam::PtrList<Foam::volScalarField> volScalarFieldCopies_;

    //- Checkpointed volVectorField fields
    Foam::PtrList<Foam::volVectorField> volVectorFields_;

    //- Checkpointed volVectorField fields (copies)
    Foam::PtrList<Foam::volVectorField> volVectorFieldCopies_;

    //- Checkpointed surfaceScalarField fields (copies)
    Foam::PtrList<Foam::surfaceScalarField> surfaceScalarFieldCopies_;

    //- Checkpointed surfaceVectorField fields (copies)
    Foam::PtrList<Foam::surfaceVectorField> surfaceVectorFieldCopies_;

    //- Checkpointed pointScalarField fields (copies)
    Foam::PtrList<Foam::pointScalarField> pointScalarFieldCopies_;

    //- Checkpointed pointVectorField fields (copies)
    Foam::PtrList<Foam::pointVectorField> pointVectorFieldCopies_;

    //- Checkpointed volTensorField fields (copies)
    Foam::PtrList<Foam::volTensorField> volTensorFieldCopies_;

    //- Checkpointed surfaceTensorField fields (copies)
    Foam::PtrList<Foam::surfaceTensorField> surfaceTensorFieldCopies_;

    //- Checkpointed pointTensorField fields (copies)
    Foam::PtrList<Foam::pointTensorField> pointTensorFieldCopies_;

    //- Checkpointed volSymmTensorField fields (copies)
    Foam::PtrList<Foam::volSymmTensorField> volSymmTensorFieldCopies_;


    // NOTE: Declare additional vectors for any other types required.

    // Configuration

    //- Read the adapter's configuration file
    bool configFileRead();

    //- Check the adapter's configuration file
    bool configFileCheck(const std::string adapterConfigFileName);

    // Methods communicating with preCICE

    //- Initialize preCICE and exchange the first data
    void initialize();

    //- Advance preCICE
    void advance();

    //- Write the coupling data at each interface
    void writeCouplingData();

    //- Determine if a checkpoint must be read
    bool isReadCheckpointRequired();

    //- Determine if a checkpoint must be written
    bool isWriteCheckpointRequired();

    //- Tell preCICE that the checkpoint has been read
    void fulfilledReadCheckpoint();

    //- Tell preCICE that the checkpoint has been written
    void fulfilledWriteCheckpoint();

    // Methods for checkpointing

    //- Configure the mesh checkpointing
    void setupMeshCheckpointing();

    //- Configure the mesh checkpointing
    void setupMeshVolCheckpointing();

    //- Configure the checkpointing
    void setupCheckpointing();

    //- Make a copy of the runTime object
    void storeCheckpointTime();

    //- Restore the copy of the runTime object
    void reloadCheckpointTime();

    //- Store the locations of the mesh points
    void storeMeshPoints();

    //- Restore the locations of the mesh points
    void reloadMeshPoints();

    // Add mesh checkpoint fields, depending on the type
    //- Add a surfaceScalarField mesh field
    void addMeshCheckpointField(surfaceScalarField& field);

    //- Add a surfaceVectorField mesh field
    void addMeshCheckpointField(surfaceVectorField& field);

    //- Add a volVectorField mesh field
    void addMeshCheckpointField(volVectorField& field);

    // TODO V0 and V00 checkpointed field.
    //- Add the V0 and V00 checkpoint fields
    void addVolCheckpointField(volScalarField::Internal& field);
    // void addVolCheckpointFieldBuffer(volScalarField::Internal & field);

    // Add checkpoint fields, depending on the type

    //- Add a volScalarField to checkpoint
    void addCheckpointField(volScalarField* field);

    //- Add a volVectorField to checkpoint
    void addCheckpointField(volVectorField* field);

    //- Add a surfaceScalarField to checkpoint
    void addCheckpointField(surfaceScalarField* field);

    //- Add a surfaceVectorField to checkpoint
    void addCheckpointField(surfaceVectorField* field);

    //- Add a pointScalarField to checkpoint
    void addCheckpointField(pointScalarField* field);

    //- Add a pointVectorField to checkpoint
    void addCheckpointField(pointVectorField* field);

    // NOTE: Add here methods to add other object types to checkpoint,
    // if needed.

    //- Add a volTensorField to checkpoint
    void addCheckpointField(volTensorField* field);

    //- Add a surfaceTensorField to checkpoint
    void addCheckpointField(surfaceTensorField* field);

    //- Add a pointTensorField to checkpoint
    void addCheckpointField(pointTensorField* field);

    //- Add a volSymmTensorField to checkpoint
    void addCheckpointField(volSymmTensorField* field);

    //- Read the checkpoint - restore the mesh fields and time
    void readMeshCheckpoint();

    //- Read the checkpoint - restore the fields and time
    void readCheckpoint();

    //- Write the checkpoint - store the mesh fields and time
    void writeMeshCheckpoint();

    //- Write the checkpoint - store the fields and time
    void writeCheckpoint();

    // TODO Probably these can be included to the mesh checkpoints.
    //- Read the volume checkpoint - restore the mesh volume fields
    void readVolCheckpoint();

    //- Write the volume checkpoint to a buffer - restore the mesh volume fields
    void writeVolCheckpoint();

    //- Store list of fields of a type
    template<class FieldType>
    void storeCheckpointFields(PtrList<FieldType>& storedFields);

    //- Load list of fields of a type
    template<class FieldType>
    void loadCheckpointFields(const PtrList<FieldType>& storedFields);

    //- Destroy the preCICE interface and delete the allocated
    //  memory in a proper way. Called by the destructor.
    void teardown();

public:
    // Methods called by the functionObject

    //- Constructor
    Adapter(const Foam::Time& runTime, const Foam::fvMesh& mesh);

    //- Setup the adapter's configuration
    void configure();

    //- Determine if the coupling timestep has been completed
    bool isCouplingTimeWindowComplete();

    //- Determine if the coupling is still happening
    bool isCouplingOngoing();

    //- Read the coupling data at each interface
    void readCouplingData();

    //- Write coupling data and advance
    void execute();

    //- Called by the functionObject's adjustTimeStep()
    scalar getMaxTimeStep();

    //- Finalize and destroy preCICE
    void finalize();

    //- Destructor
    ~Adapter();
};

}

#include "AdapterTemplates.C"

#endif
