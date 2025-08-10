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

\*---------------------------------------------------------------------------*/

#include "Adapter.H"
#include "Interface.H"
#include "Utilities.H"

#include "db/IOstreams/IOstreams.H"

using namespace Foam;

preciceAdapter::Adapter::Adapter(const Time& runTime, const fvMesh& mesh)
: runTime_(runTime),
  mesh_(mesh)
{
    adapterInfo("Loaded the FOAM-preCICE adapter", "info");

    return;
}

bool preciceAdapter::Adapter::configFileRead()
{
    SETUP_TIMER();
    adapterInfo("Reading preciceDict...", "info");

    // TODO: static is just a quick workaround to be able
    // to find the dictionary also out of scope (e.g. in KappaEffective).
    // We need a better solution.
    static IOdictionary preciceDict(
        IOobject(
            "preciceDict",
            runTime_.system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE));

    // Read and display the preCICE configuration file name
    preciceDict.lookup("preciceConfig") >> preciceConfigFilename_;
    DEBUG(adapterInfo("  precice-config-file : " + preciceConfigFilename_));

    // Read and display the participant name
    preciceDict.lookup("participant") >> participantName_;
    DEBUG(adapterInfo("  participant name    : " + participantName_));

    // Read and display the list of modules
    DEBUG(adapterInfo("  modules requested   : "));
    wordList modules;
    preciceDict.lookup("modules") >> modules;
    for (const auto& module : modules)
    {
        DEBUG(adapterInfo("  - " + module + "\n"));

        // Set the modules switches
        if (module == "CHT")
        {
            CHTenabled_ = true;
        }

        if (module == "FSI")
        {
            FSIenabled_ = true;
        }

        if (module == "FF")
        {
            FFenabled_ = true;
        }
    }

    // Every interface is a subdictionary of "interfaces",
    // each with an arbitrary name. Read all of them and create
    // a list (here: pointer) of dictionaries.
    const auto interfaceDict = preciceDict.subOrEmptyDict("interfaces");
    DEBUG(adapterInfo("  interfaces : "));

    // Check if we found any interfaces
    // and get the details of each interface
    if (interfaceDict.empty())
    {
        adapterInfo("  Empty list of interfaces", "warning");
        return false;
    }
    else
    {
        for (const entry& interfaceDictEntry : interfaceDict)
        {
            if (interfaceDictEntry.isDict())
            {
                const dictionary& interfaceDict = interfaceDictEntry.dict();
                struct InterfaceConfig interfaceConfig;

                word meshName;
                interfaceDict.lookup("mesh") >> meshName;
                interfaceConfig.meshName = meshName;
                DEBUG(adapterInfo("  - mesh         : " + interfaceConfig.meshName));

                // By default, assume "faceCenters" as locationsType
                interfaceConfig.locationsType = interfaceDict.lookupOrDefault<word>("locations", "faceCenters");
                DEBUG(adapterInfo("    locations    : " + interfaceConfig.locationsType));

                // By default, assume that no mesh connectivity is required (i.e. no nearest-projection mapping)
                interfaceConfig.meshConnectivity = interfaceDict.lookupOrDefault<bool>("connectivity", false);
                // Mesh connectivity only makes sense in case of faceNodes, check and raise a warning otherwise
                if (interfaceConfig.meshConnectivity && interfaceConfig.locationsType == "faceCenters")
                {
                    DEBUG(adapterInfo("Mesh connectivity is not supported for faceCenters. \n"
                                        "Please configure the desired interface with the locationsType faceNodes. \n"
                                        "Have a look in the adapter documentation for detailed information.",
                                        "warning"));
                    return false;
                }
                DEBUG(adapterInfo("    connectivity : " + std::to_string(interfaceConfig.meshConnectivity)));

                DEBUG(adapterInfo("    patches      : "));
                wordList patches;
                interfaceDict.lookup("patches") >> patches;
                for (auto patch : patches)
                {
                    interfaceConfig.patchNames.push_back(patch);
                    DEBUG(adapterInfo("      - " + patch));
                }

                DEBUG(adapterInfo("    writeData    : "));
                wordList writeData;
                interfaceDict.lookup("writeData") >> writeData;
                for (auto writeDatum : writeData)
                {
                    interfaceConfig.writeData.push_back(writeDatum);
                    DEBUG(adapterInfo("      - " + writeDatum));
                }

                DEBUG(adapterInfo("    readData     : "));
                wordList readData;
                interfaceDict.lookup("readData") >> readData;
                for (auto readDatum : readData)
                {
                    interfaceConfig.readData.push_back(readDatum);
                    DEBUG(adapterInfo("      - " + readDatum));
                }
                interfacesConfig_.push_back(interfaceConfig);
            }
        }
    }

    // NOTE: set the switch for your new module here

    // If the CHT module is enabled, create it, read the
    // CHT-specific options and configure it.
    if (CHTenabled_)
    {
        CHT_ = new CHT::ConjugateHeatTransfer(mesh_);
        if (!CHT_->configure(preciceDict))
        {
            return false;
        }
    }

    // If the FSI module is enabled, create it, read the
    // FSI-specific options and configure it.
    if (FSIenabled_)
    {
        // Check for unsupported FSI with meshConnectivity
        for (uint i = 0; i < interfacesConfig_.size(); i++)
        {
            if (interfacesConfig_.at(i).meshConnectivity == true)
            {
                adapterInfo(
                    "You have requested mesh connectivity (most probably for nearest-projection mapping) "
                    "and you have enabled the FSI module. "
                    "Mapping with connectivity information is not implemented for FSI, only for CHT-related fields. "
                    "warning");
                return false;
            }
        }

        FSI_ = new FSI::FluidStructureInteraction(mesh_, runTime_);
        if (!FSI_->configure(preciceDict))
        {
            return false;
        }
    }

    if (FFenabled_)
    {
        FF_ = new FF::FluidFluid(mesh_);
        if (!FF_->configure(preciceDict))
        {
            return false;
        }
    }

    // NOTE: Create your module and read any options specific to it here

    if (!CHTenabled_ && !FSIenabled_ && !FFenabled_) // NOTE: Add your new switch here
    {
        adapterInfo("No module is enabled.", "error-deferred");
        return false;
    }

    // TODO: Loading modules should be implemented in more general way,
    // in order to avoid code duplication. See issue #16 on GitHub.

    ACCUMULATE_TIMER(timeInConfigRead_);

    return true;
}

void preciceAdapter::Adapter::configure()
{
    // Read the adapter's configuration file
    if (!configFileRead())
    {
        adapterInfo(
            "There was a problem while reading the configuration file. "
            "See the log for details.",
            "error");
    }

    // Construct preCICE
    SETUP_TIMER();
    DEBUG(adapterInfo("Creating the preCICE solver interface..."));
    DEBUG(adapterInfo("  Number of processes: " + std::to_string(Pstream::nProcs())));
    DEBUG(adapterInfo("  MPI rank: " + std::to_string(Pstream::myProcNo())));
    precice_ = new precice::SolverInterface(participantName_, preciceConfigFilename_, Pstream::myProcNo(), Pstream::nProcs());
    DEBUG(adapterInfo("  preCICE solver interface was created."));

    ACCUMULATE_TIMER(timeInPreciceConstruct_);

    // Create interfaces
    REUSE_TIMER();
    DEBUG(adapterInfo("Creating interfaces..."));
    for (uint i = 0; i < interfacesConfig_.size(); i++)
    {
        std::string namePointDisplacement = FSIenabled_ ? FSI_->getPointDisplacementFieldName() : "default";
        std::string nameCellDisplacement = FSIenabled_ ? FSI_->getCellDisplacementFieldName() : "default";
        bool restartFromDeformed = FSIenabled_ ? FSI_->isRestartingFromDeformed() : false;

        Interface* interface = new Interface(*precice_, mesh_, interfacesConfig_.at(i).meshName, interfacesConfig_.at(i).locationsType, interfacesConfig_.at(i).patchNames, interfacesConfig_.at(i).meshConnectivity, restartFromDeformed, namePointDisplacement, nameCellDisplacement);
        interfaces_.push_back(interface);
        DEBUG(adapterInfo("Interface created on mesh " + interfacesConfig_.at(i).meshName));

        DEBUG(adapterInfo("Adding coupling data writers..."));
        for (uint j = 0; j < interfacesConfig_.at(i).writeData.size(); j++)
        {
            std::string dataName = interfacesConfig_.at(i).writeData.at(j);

            unsigned int inModules = 0;

            // Add CHT-related coupling data writers
            if (CHTenabled_ && CHT_->addWriters(dataName, interface))
            {
                inModules++;
            }

            // Add FSI-related coupling data writers
            if (FSIenabled_ && FSI_->addWriters(dataName, interface))
            {
                inModules++;
            }

            // Add FF-related coupling data writers
            if (FFenabled_ && FF_->addWriters(dataName, interface))
            {
                inModules++;
            }

            if (inModules == 0)
            {
                adapterInfo("I don't know how to write \"" + dataName
                                + "\". Maybe this is a typo or maybe you need to enable some adapter module?",
                            "error-deferred");
            }
            else if (inModules > 1)
            {
                adapterInfo("It looks like more than one modules can write \"" + dataName
                                + "\" and I don't know how to choose. Try disabling one of the modules.",
                            "error-deferred");
            }

            // NOTE: Add any coupling data writers for your module here.
        } // end add coupling data writers

        DEBUG(adapterInfo("Adding coupling data readers..."));
        for (uint j = 0; j < interfacesConfig_.at(i).readData.size(); j++)
        {
            std::string dataName = interfacesConfig_.at(i).readData.at(j);

            unsigned int inModules = 0;

            // Add CHT-related coupling data readers
            if (CHTenabled_ && CHT_->addReaders(dataName, interface)) inModules++;

            // Add FSI-related coupling data readers
            if (FSIenabled_ && FSI_->addReaders(dataName, interface)) inModules++;

            // Add FF-related coupling data readers
            if (FFenabled_ && FF_->addReaders(dataName, interface)) inModules++;

            if (inModules == 0)
            {
                adapterInfo("I don't know how to read \"" + dataName
                                + "\". Maybe this is a typo or maybe you need to enable some adapter module?",
                            "error-deferred");
            }
            else if (inModules > 1)
            {
                adapterInfo("It looks like more than one modules can read \"" + dataName
                                + "\" and I don't know how to choose. Try disabling one of the modules.",
                            "error-deferred");
            }

            // NOTE: Add any coupling data readers for your module here.
        } // end add coupling data readers

        // Create the interface's data buffer
        interface->createBuffer();
    }
    ACCUMULATE_TIMER(timeInMeshSetup_);

    // Initialize preCICE and exchange the first coupling data
    initialize();

    // Read the received coupling data
    readCouplingData();

    // If checkpointing is required, specify the checkpointed fields
    // and write the first checkpoint
    if (isWriteCheckpointRequired())
    {
        checkpointing_ = true;

        // Setup the checkpointing (find and add fields to checkpoint)
        setupCheckpointing();

        // Write checkpoint (for the first iteration)
        writeCheckpoint();
        fulfilledWriteCheckpoint();
    }
}

void preciceAdapter::Adapter::execute()
{
    // The solver has already solved the equations for this timestep.
    // Now call the adapter's methods to perform the coupling.

    // TODO add a function which checks if all fields are checkpointed.
    // if (ncheckpointed is nregisterdobjects.)

    // Write the coupling data in the buffer
    writeCouplingData();

    // Advance preCICE
    advance();

    // Read checkpoint if required
    if (isReadCheckpointRequired())
    {
        readCheckpoint();
        fulfilledReadCheckpoint();
    }

    // Write checkpoint if required
    if (isWriteCheckpointRequired())
    {
        writeCheckpoint();
        fulfilledWriteCheckpoint();
    }

    if (!isCouplingOngoing())
    {
        adapterInfo("The coupling completed.", "info");
        finalize();
    }
}

void preciceAdapter::Adapter::readCouplingData()
{
    SETUP_TIMER();
    DEBUG(adapterInfo("Reading coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->readCouplingData();
    }

    ACCUMULATE_TIMER(timeInRead_);

    return;
}

void preciceAdapter::Adapter::writeCouplingData()
{
    SETUP_TIMER();
    DEBUG(adapterInfo("Writing coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->writeCouplingData();
    }

    ACCUMULATE_TIMER(timeInWrite_);

    return;
}

void preciceAdapter::Adapter::initialize()
{
    DEBUG(adapterInfo("Initializing the preCICE solver interface..."));
    SETUP_TIMER();
    timestepPrecice_ = precice_->initialize();
    ACCUMULATE_TIMER(timeInInitialize_);

    preciceInitialized_ = true;

    if (precice_->isActionRequired(precice::constants::actionWriteInitialData()))
    {
        writeCouplingData();
        precice_->markActionFulfilled(precice::constants::actionWriteInitialData());
    }

    DEBUG(adapterInfo("Initializing preCICE data..."));
    REUSE_TIMER();
    precice_->initializeData();
    ACCUMULATE_TIMER(timeInInitializeData_);

    adapterInfo("preCICE was configured and initialized", "info");

    return;
}

void preciceAdapter::Adapter::finalize()
{
    if (NULL != precice_ && preciceInitialized_ && !isCouplingOngoing())
    {
        DEBUG(adapterInfo("Finalizing the preCICE solver interface..."));

        // Finalize the preCICE solver interface
        SETUP_TIMER();
        precice_->finalize();
        ACCUMULATE_TIMER(timeInFinalize_);

        preciceInitialized_ = false;
    }

    // Delete the solver interface and all the related data
    teardown();
}

void preciceAdapter::Adapter::advance()
{
    DEBUG(adapterInfo("Advancing preCICE..."));

    SETUP_TIMER();
    timestepPrecice_ = precice_->advance(runTime_.deltaT().value());
    ACCUMULATE_TIMER(timeInAdvance_);

    return;
}

scalar preciceAdapter::Adapter::getMaxTimeStep()
{
    DEBUG(adapterInfo("Adjusting the solver's timestep..."));

    return timestepPrecice_;
}

bool preciceAdapter::Adapter::isCouplingOngoing()
{
    bool isCouplingOngoing = false;

    // If the coupling ends before the solver,
    // the solver would try to access this method again,
    // giving a segmentation fault if precice_
    // was not available.
    if (NULL != precice_)
    {
        isCouplingOngoing = precice_->isCouplingOngoing();
    }

    return isCouplingOngoing;
}

bool preciceAdapter::Adapter::isCouplingTimeWindowComplete()
{
    return precice_->isTimeWindowComplete();
}

bool preciceAdapter::Adapter::isReadCheckpointRequired()
{
    return precice_->isActionRequired(precice::constants::actionReadIterationCheckpoint());
}

bool preciceAdapter::Adapter::isWriteCheckpointRequired()
{
    return precice_->isActionRequired(precice::constants::actionWriteIterationCheckpoint());
}

void preciceAdapter::Adapter::fulfilledReadCheckpoint()
{
    precice_->markActionFulfilled(precice::constants::actionReadIterationCheckpoint());
}

void preciceAdapter::Adapter::fulfilledWriteCheckpoint()
{
    precice_->markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
}

void preciceAdapter::Adapter::storeCheckpointTime()
{
    couplingIterationTimeIndex_ = runTime_.timeIndex();
    couplingIterationTimeValue_ = runTime_.value();
    DEBUG(adapterInfo("Stored time value t = " + std::to_string(runTime_.value())));
}

void preciceAdapter::Adapter::reloadCheckpointTime()
{
    const_cast<Time&>(runTime_).setTime(couplingIterationTimeValue_, couplingIterationTimeIndex_);
    // TODO also reset the current iteration?!
    DEBUG(adapterInfo("Reloaded time value t = " + std::to_string(runTime_.value())));
}

void preciceAdapter::Adapter::storeMeshPoints()
{
    DEBUG(adapterInfo("Storing mesh points..."));
    meshPoints_ = mesh_.points();
    oldMeshPoints_ = mesh_.oldPoints();

    /*
    // TODO  This is only required for subcycling. It should not be called when not subcycling!!
    // Add a bool 'subcycling' which can be evaluated every timestep.
    if (!oldVolsStored && mesh_.foundObject<volScalarField::Internal>("V00")) // For Ddt schemes which use one previous timestep
    {
        setupMeshVolCheckpointing();
        oldVolsStored = true;
    }
    // Update any volume fields from the buffer to the checkpointed values (if already exists.)
    */

    DEBUG(adapterInfo("Stored mesh points."));
    if (mesh_.moving())
    {
        if (!meshCheckPointed)
        {
            // Set up the checkpoint for the mesh flux: meshPhi
            setupMeshCheckpointing();
            meshCheckPointed = true;
        }
        writeMeshCheckpoint();
        writeVolCheckpoint(); // Does not write anything unless subcycling.
    }
}

void preciceAdapter::Adapter::reloadMeshPoints()
{
    if (!mesh_.moving())
    {
        DEBUG(adapterInfo("Mesh points not moved as the mesh is not moving"));
        return;
    }

    // In Foam::polyMesh::movePoints.
    // TODO: The function movePoints overwrites the pointer to the old mesh.
    // Therefore, if you revert the mesh, the oldpointer will be set to the points, which are the new values.
    DEBUG(adapterInfo("Moving mesh points to their previous locations..."));

    // TODO
    // Switch oldpoints on for pure physics. (is this required?). Switch off for better mesh deformation capabilities?
    // const_cast<pointField&>(mesh_.points()) = oldMeshPoints_;
    const_cast<fvMesh&>(mesh_).movePoints(meshPoints_);

    DEBUG(adapterInfo("Moved mesh points to their previous locations."));

    // TODO The if statement can be removed in this case, but it is still included for clarity
    if (meshCheckPointed)
    {
        readMeshCheckpoint();
    }

    /*  // TODO This part should only be used when sybcycling. See the description in 'storeMeshPoints()'
        // The if statement can be removed in this case, but it is still included for clarity
    if (oldVolsStored)
    {
        readVolCheckpoint();
    }
    */
}

void preciceAdapter::Adapter::setupMeshCheckpointing()
{
    // The other mesh <type>Fields:
    //      C
    //      Cf
    //      Sf
    //      magSf
    //      delta
    // are updated by the function fvMesh::movePoints. Only the meshPhi needs checkpointing.
    DEBUG(adapterInfo("Creating a list of the mesh checkpointed fields..."));

    // Add meshPhi to the checkpointed fields
    addMeshCheckpointField(
        const_cast<surfaceScalarField&>(
            mesh_.phi()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.phi().name() + " in the list of checkpointed fields.");
#endif
}

void preciceAdapter::Adapter::setupMeshVolCheckpointing()
{
    DEBUG(adapterInfo("Creating a list of the mesh volume checkpointed fields..."));
    // Add the V0 and the V00 to the list of checkpointed fields.
    // For V0
    addVolCheckpointField(
        const_cast<volScalarField::Internal&>(
            mesh_.V0()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V0().name() + " in the list of checkpointed fields.");
#endif
    // For V00
    addVolCheckpointField(
        const_cast<volScalarField::Internal&>(
            mesh_.V00()));
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V00().name() + " in the list of checkpointed fields.");
#endif

    // Also add the buffer fields.
    // TODO For V0
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V0()
        )
    ); */
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V0().name() + " in the list of buffer checkpointed fields.");
#endif
    // TODO For V00
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V00()
        )
    );*/
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V00().name() + " in the list of buffer checkpointed fields.");
#endif
}


void preciceAdapter::Adapter::setupCheckpointing()
{
    SETUP_TIMER();

    // Add fields in the checkpointing list - sorted for parallel consistency
    DEBUG(adapterInfo("Adding in checkpointed fields..."));

#undef doLocalCode
#define doLocalCode(GeomField)                                           \
    /* Checkpoint registered GeomField objects */                        \
    for (const word& obj : mesh_.sortedNames<GeomField>())               \
    {                                                                    \
        addCheckpointField(mesh_.thisDb().getObjectPtr<GeomField>(obj)); \
        DEBUG(adapterInfo("Checkpoint " + obj + " : " #GeomField));      \
    }

    doLocalCode(volScalarField);
    doLocalCode(volVectorField);
    doLocalCode(volTensorField);
    doLocalCode(volSymmTensorField);

    doLocalCode(surfaceScalarField);
    doLocalCode(surfaceVectorField);
    doLocalCode(surfaceTensorField);

    doLocalCode(pointScalarField);
    doLocalCode(pointVectorField);
    doLocalCode(pointTensorField);

    // NOTE: Add here other object types to checkpoint, if needed.

#undef doLocalCode

    ACCUMULATE_TIMER(timeInCheckpointingSetup_);
}


// All mesh checkpointed fields

void preciceAdapter::Adapter::addMeshCheckpointField(surfaceScalarField& field)
{
    meshSurfaceScalarFieldCopies_.append(new surfaceScalarField(field));
}

void preciceAdapter::Adapter::addMeshCheckpointField(surfaceVectorField& field)
{
    meshSurfaceVectorFieldCopies_.append(new surfaceVectorField(field));
}

void preciceAdapter::Adapter::addMeshCheckpointField(volVectorField& field)
{
    meshVolVectorFieldCopies_.append(new volVectorField(field));
}

// TODO Internal field for the V0 (volume old) and V00 (volume old-old) fields
void preciceAdapter::Adapter::addVolCheckpointField(volScalarField::Internal& field)
{
    volScalarInternalFieldCopies_.append(new volScalarField::Internal(field));
}


void preciceAdapter::Adapter::addCheckpointField(volScalarField* field)
{
    if (field)
    {
        volScalarFieldCopies_.append(new volScalarField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(volVectorField* field)
{
    if (field)
    {
        volVectorFieldCopies_.append(new volVectorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(surfaceScalarField* field)
{
    if (field)
    {
        surfaceScalarFieldCopies_.append(new surfaceScalarField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(surfaceVectorField* field)
{
    if (field)
    {
        surfaceVectorFieldCopies_.append(new surfaceVectorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(pointScalarField* field)
{
    if (field)
    {
        pointScalarFieldCopies_.append(new pointScalarField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(pointVectorField* field)
{
    if (field)
    {
        pointVectorFieldCopies_.append(new pointVectorField(*field));
        // TODO: Old time
        // pointVectorFieldCopiesOld_.append(new pointVectorField(field->oldTime()));
    }
}

void preciceAdapter::Adapter::addCheckpointField(volTensorField* field)
{
    if (field)
    {
        volTensorFieldCopies_.append(new volTensorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(surfaceTensorField* field)
{
    if (field)
    {
        surfaceTensorFieldCopies_.append(new surfaceTensorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(pointTensorField* field)
{
    if (field)
    {
        pointTensorFieldCopies_.append(new pointTensorField(*field));
    }
}

void preciceAdapter::Adapter::addCheckpointField(volSymmTensorField* field)
{
    if (field)
    {
        volSymmTensorFieldCopies_.append(new volSymmTensorField(*field));
    }
}


// NOTE: Add here methods to add other object types to checkpoint, if needed.

void preciceAdapter::Adapter::readCheckpoint()
{
    SETUP_TIMER();

    // TODO: To increase efficiency: only the oldTime() fields of the quantities which are used in the time
    //  derivative are necessary. (In general this is only the velocity). Also old information of the mesh
    //  is required.
    //  Therefore, loading the oldTime() and oldTime().oldTime() fields for the other fields can be excluded
    //  for efficiency.
    DEBUG(adapterInfo("Reading a checkpoint..."));

    // Reload the meshPoints (if FSI is enabled)
    if (FSIenabled_)
    {
        reloadMeshPoints();
    }

    loadCheckpointFields(volScalarFieldCopies_);
    loadCheckpointFields(volVectorFieldCopies_);
    loadCheckpointFields(surfaceScalarFieldCopies_);
    loadCheckpointFields(surfaceVectorFieldCopies_);
    loadCheckpointFields(pointScalarFieldCopies_);
    loadCheckpointFields(pointVectorFieldCopies_);
    // TODO Evaluate if all the tensor fields need to be in here.
    loadCheckpointFields(volTensorFieldCopies_);
    loadCheckpointFields(surfaceTensorFieldCopies_);
    loadCheckpointFields(pointTensorFieldCopies_);
    loadCheckpointFields(volSymmTensorFieldCopies_);

    // NOTE: Add here other field types to read, if needed.

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Checkpoint was read. Time = " + std::to_string(runTime_.value()));
#endif

    ACCUMULATE_TIMER(timeInCheckpointingRead_);
}


void preciceAdapter::Adapter::writeCheckpoint()
{
    SETUP_TIMER();

    DEBUG(adapterInfo("Writing a checkpoint..."));

    // Store the runTime
    storeCheckpointTime();

    // Store the meshPoints (if FSI is enabled)
    if (FSIenabled_)
    {
        storeMeshPoints();
    }

    storeCheckpointFields(volScalarFieldCopies_);
    storeCheckpointFields(volVectorFieldCopies_);
    storeCheckpointFields(volTensorFieldCopies_);
    storeCheckpointFields(volSymmTensorFieldCopies_);
    storeCheckpointFields(surfaceScalarFieldCopies_);
    storeCheckpointFields(surfaceVectorFieldCopies_);
    storeCheckpointFields(surfaceTensorFieldCopies_);
    storeCheckpointFields(pointScalarFieldCopies_);
    storeCheckpointFields(pointVectorFieldCopies_);
    storeCheckpointFields(pointTensorFieldCopies_);

    // NOTE: Add here other types to write, if needed.

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif

    ACCUMULATE_TIMER(timeInCheckpointingWrite_);
}

void preciceAdapter::Adapter::readMeshCheckpoint()
{
    DEBUG(adapterInfo("Reading a mesh checkpoint..."));

    // TODO only the meshPhi field is here, which is a surfaceScalarField. The other fields can be removed.
    //  Reload all the fields of type mesh surfaceScalarField
    loadCheckpointFields(meshSurfaceScalarFieldCopies_);
    loadCheckpointFields(meshSurfaceVectorFieldCopies_);
    loadCheckpointFields(meshVolVectorFieldCopies_);

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh checkpoint was read. Time = " + std::to_string(runTime_.value()));
#endif
}

void preciceAdapter::Adapter::writeMeshCheckpoint()
{
    DEBUG(adapterInfo("Writing a mesh checkpoint..."));

    loadCheckpointFields(meshSurfaceScalarFieldCopies_);
    loadCheckpointFields(meshSurfaceVectorFieldCopies_);
    loadCheckpointFields(meshVolVectorFieldCopies_);

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif
}

// TODO for the volumes of the mesh, check this part for subcycling.
void preciceAdapter::Adapter::readVolCheckpoint()
{
    DEBUG(adapterInfo("Reading the mesh volumes checkpoint..."));

    // Reload all the fields of type mesh volVectorField::Internal
    for (const volScalarField::Internal& fc : volScalarInternalFieldCopies_)
    {
        // Load the volume field
        volScalarField::Internal& f =
            fc.db().template lookupObjectRef<volScalarField::Internal>
            (
                fc.name()
            );
        f = fc;
        // There are no old times for the internal fields.
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh volumes were read. Time = " + std::to_string(runTime_.value()));
#endif
}

void preciceAdapter::Adapter::writeVolCheckpoint()
{
    DEBUG(adapterInfo("Writing a mesh volumes checkpoint..."));

    // Store all the fields of type mesh volScalarField::Internal
    for (volScalarField::Internal& f : volScalarInternalFieldCopies_)
    {
        f = f.db().template lookupObject<volScalarField::Internal>(f.name());
    }

#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Mesh volumes checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
#endif
}


void preciceAdapter::Adapter::teardown()
{
    // If the solver interface was not deleted before, delete it now.
    // Normally it should be deleted when isCouplingOngoing() becomes false.
    if (NULL != precice_)
    {
        DEBUG(adapterInfo("Destroying the preCICE solver interface..."));
        delete precice_;
        precice_ = NULL;
    }

    // Delete the preCICE solver interfaces
    if (interfaces_.size() > 0)
    {
        DEBUG(adapterInfo("Deleting the interfaces..."));
        for (uint i = 0; i < interfaces_.size(); i++)
        {
            delete interfaces_.at(i);
        }
        interfaces_.clear();
    }

    // Delete the copied fields for checkpointing
    if (checkpointing_)
    {
        DEBUG(adapterInfo("Deleting the checkpoints... "));

        // Fields
        // volScalarFields
        volScalarFieldCopies_.clear();

        // volVector
        volVectorFieldCopies_.clear();

        // surfaceScalar
        surfaceScalarFieldCopies_.clear();

        // surfaceVector
        surfaceVectorFieldCopies_.clear();

        // pointScalar
        pointScalarFieldCopies_.clear();

        // pointVector
        pointVectorFieldCopies_.clear();

        // Mesh fields
        // meshSurfaceScalar
        meshSurfaceScalarFieldCopies_.clear();

        // meshSurfaceVector
        meshSurfaceVectorFieldCopies_.clear();

        // meshVolVector
        meshVolVectorFieldCopies_.clear();

        // TODO for the internal volume
        //  volScalarInternal
        volScalarInternalFieldCopies_.clear();

        // volTensorField
        volTensorFieldCopies_.clear();

        // surfaceTensorField
        surfaceTensorFieldCopies_.clear();

        // pointTensorField
        pointTensorFieldCopies_.clear();

        // volSymmTensor
        volSymmTensorFieldCopies_.clear();

        // NOTE: Add here delete for other types, if needed

        checkpointing_ = false;
    }

    // Delete the CHT module
    if (NULL != CHT_)
    {
        DEBUG(adapterInfo("Destroying the CHT module..."));
        delete CHT_;
        CHT_ = NULL;
    }

    // Delete the FSI module
    if (NULL != FSI_)
    {
        DEBUG(adapterInfo("Destroying the FSI module..."));
        delete FSI_;
        FSI_ = NULL;
    }

    // Delete the FF module
    if (NULL != FF_)
    {
        DEBUG(adapterInfo("Destroying the FF module..."));
        delete FF_;
        FF_ = NULL;
    }

    // NOTE: Delete your new module here
}

preciceAdapter::Adapter::~Adapter()
{
    // Delete data if not already done
    teardown();

    adapterInfo("The FOAM preCICE adapter is based on the OpenFOAM-preCICE adapter. "
                "Please see https://precice.org/adapter-openfoam-overview.html.",
                "info");

    TIMING_MODE(
        // Continuing the output started in the destructor of preciceAdapterFunctionObject
        Info<< "Time exclusively in the adapter: " << (timeInConfigRead_ + timeInMeshSetup_ + timeInCheckpointingSetup_ + timeInWrite_ + timeInRead_ + timeInCheckpointingWrite_ + timeInCheckpointingRead_).str() << nl;
        Info<< "  (S) reading preciceDict:       " << timeInConfigRead_.str() << nl;
        Info<< "  (S) constructing preCICE:      " << timeInPreciceConstruct_.str() << nl;
        Info<< "  (S) setting up the interfaces: " << timeInMeshSetup_.str() << nl;
        Info<< "  (S) setting up checkpointing:  " << timeInCheckpointingSetup_.str() << nl;
        Info<< "  (I) writing data:              " << timeInWrite_.str() << nl;
        Info<< "  (I) reading data:              " << timeInRead_.str() << nl;
        Info<< "  (I) writing checkpoints:       " << timeInCheckpointingWrite_.str() << nl;
        Info<< "  (I) reading checkpoints:       " << timeInCheckpointingRead_.str() << nl;
        Info<< "  (I) writing OpenFOAM results:  " << timeInWriteResults_.str() << " (at the end of converged time windows)" << nl << nl;
        Info<< "Time exclusively in preCICE:     " << (timeInInitialize_ + timeInInitializeData_ + timeInAdvance_ + timeInFinalize_).str() << nl;
        Info<< "  (S) initialize():              " << timeInInitialize_.str() << nl;
        Info<< "  (S) initializeData():          " << timeInInitializeData_.str() << nl;
        Info<< "  (I) advance():                 " << timeInAdvance_.str() << nl;
        Info<< "  (I) finalize():                " << timeInFinalize_.str() << nl;
        Info<< "  These times include time waiting for other participants." << nl;
        Info<< "  See also precice-<participant>-events-summary.log." << nl;
        Info<< "-------------------------------------------------------------------------------------" << nl;)
}
