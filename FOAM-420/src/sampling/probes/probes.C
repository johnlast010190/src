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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "probes/probes.H"
#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(probes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        probes,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::probes::findElements(const fvMesh& mesh)
{
    DebugInformation<< "probes: resetting sample locations" << endl;

    //check for conflicting

    elementList_.clear();
    elementList_.setSize(size());

    faceList_.clear();
    faceList_.setSize(size());

    processor_.setSize(size());
    processor_ = -1;

    if (nearestElement_)
    {
        elementList_ = -1;
        List<List<scalar>> cellDist(Pstream::nProcs());

        cellDist[Pstream::myProcNo()].setSize(this->size());
        cellDist[Pstream::myProcNo()] = GREAT;

        List<scalar>& myCellDist = cellDist[Pstream::myProcNo()];

        forAll(*this, probei)
        {
            const vector& location = operator[](probei);
            elementList_[probei] = mesh.findNearestCell(location);
            if (elementList_[probei] != -1)
            {
                myCellDist[probei]
                    = mag(mesh.cellCentres()[elementList_[probei]] - location);
            }
        }

        Pstream::allGatherList(cellDist);

        forAll(*this, probei)
        {
            label procWithClosestCell = -1;
            scalar closestCellDist = GREAT;

            forAll(cellDist, proci)
            {
                if (cellDist[proci][probei] < closestCellDist)
                {
                    closestCellDist = cellDist[proci][probei];
                    procWithClosestCell = proci;
                }
                else if
                (
                    cellDist[proci][probei] == closestCellDist
                    && proci < procWithClosestCell
                )
                {
                    closestCellDist = cellDist[proci][probei];
                    procWithClosestCell = proci;
                }
            }

            if (Pstream::myProcNo() !=  procWithClosestCell)
            {
                elementList_[probei] = -1;
            }

            if (debug && elementList_[probei] != -1)
            {
                Pout<< "volProbes : " << this->operator[](probei)
                    << " closest cell " << elementList_[probei] << endl;

            }
        }
    }
    else
    {
        forAll(*this, probei)
        {
            const vector& location = operator[](probei);

            const label celli = mesh.findCell(location);

            elementList_[probei] = celli;
        }
    }

    forAll(*this, probei)
    {
        const label celli = elementList_[probei];
        const vector& location = operator[](probei);

        if (celli != -1)
        {
            const labelList& cellFaces = mesh.cells()[celli];

            scalar minDistance = GREAT;
            label minFaceID = -1;

            forAll(cellFaces, i)
            {
                label facei = cellFaces[i];
                vector dist = mesh.faceCentres()[facei] - location;
                if (mag(dist) < minDistance)
                {
                    minDistance = mag(dist);
                    minFaceID = facei;
                }
            }

            faceList_[probei] = minFaceID;
        }
        else
        {
            faceList_[probei] = -1;
        }

        if (debug && (elementList_[probei] != -1 || faceList_[probei] != -1))
        {
            Pout<< "probes : found point " << location
                << " in cell " << elementList_[probei]
                << " and face " << faceList_[probei] << endl;
        }
    }

    // Check if all probes have been found.
    forAll(elementList_, probei)
    {
        const vector& location = operator[](probei);
        label celli = elementList_[probei];
        label facei = faceList_[probei];

        processor_[probei] = (celli != -1 ? Pstream::myProcNo() : -1);

        // Check at least one processor with cell.
        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());
        reduce(processor_[probei], maxOp<label>());

        if (celli == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (facei == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (elementList_[probei] != -1 && elementList_[probei] != celli)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << elementList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and cell " << celli << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceList_[probei] != -1 && faceList_[probei] != facei)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and face " << facei << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}


Foam::label Foam::probes::prepare()
{
    const label nFields = classifyFields();
    wordHashSet allFieldNames = getAllFieldNames();
    const fileName outputDir = getOutputDir();
    if (newFileFormat_)
    {
        initialiseOutputFiles(allFieldNames, outputDir);
    }
    else
    {
        initialiseLegacyOutputFiles(allFieldNames, outputDir);
    }

    return nFields;
}


Foam::wordHashSet Foam::probes::getAllFieldNames() const
{
    wordHashSet allFieldNames;

    allFieldNames.insert(scalarFields_);
    allFieldNames.insert(vectorFields_);
    allFieldNames.insert(sphericalTensorFields_);
    allFieldNames.insert(symmTensorFields_);
    allFieldNames.insert(tensorFields_);

    allFieldNames.insert(surfaceScalarFields_);
    allFieldNames.insert(surfaceVectorFields_);
    allFieldNames.insert(surfaceSphericalTensorFields_);
    allFieldNames.insert(surfaceSymmTensorFields_);
    allFieldNames.insert(surfaceTensorFields_);

    DebugInformation
        << "Probing fields: " << allFieldNames << nl
        << "Probing locations: " << *this << nl
        << endl;

    return allFieldNames;
}


const Foam::fileName Foam::probes::getOutputDir() const
{
    fileName probeDir;
    fileName probeSubDir = name();

    if (mesh_.name() != polyMesh::defaultRegion)
    {
        probeSubDir = mesh_.name()/probeSubDir;
    }
    probeSubDir = "postProcessing"/probeSubDir/mesh_.time().timeName();

    if (Pstream::parRun())
    {
        // Put in undecomposed case
        // (Note: gives problems for distributed data running)
        probeDir = mesh_.time().path()/".."/probeSubDir;
    }
    else
    {
        probeDir = mesh_.time().path()/probeSubDir;
    }
    // Remove ".."
    probeDir.clean();

    return probeDir;
}


void Foam::probes::initialiseLegacyOutputFiles
(
    wordHashSet& currentFields,
    const fileName& probeDir
)
{
    // (legacy) Delete all file handles which have keys which don't exist in the
    // currentFields HashSet
    forAllIter(HashPtrTable<OFstream>, probeFilePtrs_, iter)
    {
        if (!currentFields.erase(iter.key()))
        {
            DebugInformation<< "close probe stream: " << iter()->name() << endl;

            delete probeFilePtrs_.remove(iter);
        }
    }

    // (legacy) currentFields now just has the new fields - open streams for
    // them
    forAllConstIter(wordHashSet, currentFields, iter)
    {
        const word& fieldName = iter.key();

        // Create directory if does not exist.
        mkDir(probeDir);
        OFstream* fPtr = new OFstream(probeDir/fieldName);
        OFstream& fout = *fPtr;
        DebugInformation<< "open probe stream: " << fout.name() << endl;
        probeFilePtrs_.insert(fieldName, fPtr);
        unsigned int w = IOstream::defaultPrecision() + 7;
        forAll(*this, probei)
        {
            fout<< "# Probe " << probei << ' ' << operator[](probei);

            if (processor_[probei] == -1)
            {
                fout<< "  # Not Found";
            }
            fout<< endl;
        }
        fout<< '#' << setw(IOstream::defaultPrecision() + 6)
            << "Probe";

        forAll(*this, probei)
        {
            if (includeOutOfBounds_ || processor_[probei] != -1)
            {
                fout<< ' ' << setw(w) << probei;
                if (coorFramePtr_)
                {
                    fout<< " (location value)";
                }
            }
        }
        fout<< endl;

        fout<< '#' << setw(IOstream::defaultPrecision() + 6)
            << "Time" << endl;
    }
}

void Foam::probes::initialiseOutputFiles
(
    wordHashSet& allFieldNames,
    const fileName& probeDir
)
{
    // Delete all file writers which have keys which don't exist in allFieldNames
    forAllIter(HashPtrTable<columnatedFileWriter>, fileWriters_, iter)
    {
        if (!allFieldNames.erase(iter.key()))
        {
            DebugInformation << "Closing file writer for " << iter.key() << " probe" << endl;

            delete fileWriters_.remove(iter);
        }
    }

    for (const std::string& fieldName : allFieldNames)
    {
        // Horrible, but the only sensible solution for now - see note in the
        // header file
        columnatedFileWriter* fileWriter =
            new columnatedFileWriter
            (
                mesh_.time(),
                probeDir,
                fieldName,
                dict_
            );
        fileWriter->writeCommented("Time");
        forAll(*this, probei)
        {
            std::stringstream probeId;
            probeId << "Probe " << probei;
            if (processor_[probei] == -1)
            {
                probeId<< "  # Not Found";
            }
            else if (coorFramePtr_)
            {
                probeId << " (location value)";
            }
            fileWriter->writeDelimited(probeId.str());
        }
        fileWriter->endLine();
        fileWriters_.insert(fieldName, fileWriter);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::probes::probes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    regionFunctionObject(name, runTime, dict),
    pointField(0),
    dict_(dict),
    newFileFormat_(dict.found("outputFileFormat")),
    mesh_(runTime.lookupObject<fvMesh>(meshObr_.name())),
    loadFromFiles_(loadFromFiles),
    meshUpdated_(false),
    coorFramePtr_(nullptr),
    localOnOutput_(false),
    fieldSelection_(),
    fixedLocations_(true),
    interpolationScheme_("cell"),
    includeOutOfBounds_(true),
    nearestElement_(true)
{
    if (readFields)
    {
        read(dict);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::probes::~probes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::probes::read(const dictionary& dict)
{
    dict.lookup("probeLocations") >> *this;
    dict.lookup("fields") >> fieldSelection_;

    if (dict.found("referenceFrame"))
    {
        coorFramePtr_ = coordinateFrame::lookupNew(mesh_, dict);
        if (dict.lookupOrDefault<Switch>("definedInFrame", false))
        {
            pointField& samplePoints = *this;
            samplePoints =
                coorFramePtr_->coorSys().globalPosition(samplePoints).ref();
        }
        localOnOutput_ = dict.lookupOrDefault<Switch>("localOnOutput", false);

        // Try to read point values when simulation is re-started
        if (time_.value() != 0.0 && coorFramePtr_->isIncrementalMotion())
        {
            IOobject header
            (
                "functions" + name() + "Locations",
                obr_.parent().time().timeName(),
                "uniform",
                obr_.parent(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            );
            if (header.typeHeaderOk<IOdictionary>(true))
            {
                IOdictionary oldPointsDict(header);
                dictionary(oldPointsDict).lookup("probeLocations") >> *this;
            }
        }
    }

    dict.readIfPresent("fixedLocations", fixedLocations_);
    if (dict.readIfPresent("interpolationScheme", interpolationScheme_))
    {
        if (!fixedLocations_ && interpolationScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations.  InterpolationScheme "
                << "entry will be ignored"
                << endl;
        }
    }
    dict.readIfPresent("includeOutOfBounds", includeOutOfBounds_);
    dict.readIfPresent("nearestElement", nearestElement_);

    // Initialise cells to sample from supplied locations
    findElements(mesh_);

    prepare();

    return true;
}


bool Foam::probes::execute()
{
    //If using fixedLocations update if mesh has been updated
    if (meshUpdated_ || coorFramePtr_)
    {
        if (coorFramePtr_)
        {
            coorFramePtr_->updateState();
            pointField& samplePoints = *this;
            if (!coorFramePtr_->isIncrementalMotion())
            {
                if (time0Points_.empty())
                {
                    time0Points_ = new vectorField(*this);
                }
                else
                {
                    samplePoints = pointField(time0Points_.ref());
                }
            }
            samplePoints =
                    transformPoints
                    (
                        coorFramePtr_->transformation(),
                        samplePoints
                    );
        }
        findElements(mesh_);
        meshUpdated_ = false;
    }

    if (!empty() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


bool Foam::probes::write()
{
    if
    (
        time_.writeTime()
     && coorFramePtr_
     && coorFramePtr_->isIncrementalMotion()
    )
    {
        IOdictionary dict
        (
            IOobject
            (
                "functions" + name() + "Locations",
                obr_.parent().time().timeName(),
                "uniform",
                obr_.parent(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        dict.set("probeLocations", *this);
        dict.regIOobject::write();
    }

    return true;
}


void Foam::probes::updateMesh(const mapPolyMesh& mpm)
{
    DebugInformation<< "probes: updateMesh" << endl;

    if (&mpm.mesh() != &mesh_)
    {
        return;
    }

    if (fixedLocations_)
    {
        meshUpdated_ = true;
    }
    else
    {
        DebugInformation<< "probes: remapping sample locations" << endl;

        // 1. Update cells
        {
            DynamicList<label> elems(elementList_.size());

            const labelList& reverseMap = mpm.reverseCellMap();
            forAll(elementList_, i)
            {
                label celli = elementList_[i];
                if (celli != -1)
                {
                    label newCelli = reverseMap[celli];
                    if (newCelli == -1)
                    {
                        // cell removed
                    }
                    else if (newCelli < -1)
                    {
                        // cell merged
                        elems.append(-newCelli - 2);
                    }
                    else
                    {
                        // valid new cell
                        elems.append(newCelli);
                    }
                }
                else
                {
                    // Keep -1 elements so the size stays the same
                    elems.append(-1);
                }
            }

            elementList_.transfer(elems);
        }

        // 2. Update faces
        {
            DynamicList<label> elems(faceList_.size());

            const labelList& reverseMap = mpm.reverseFaceMap();
            forAll(faceList_, i)
            {
                label facei = faceList_[i];
                if (facei != -1)
                {
                    label newFacei = reverseMap[facei];
                    if (newFacei == -1)
                    {
                        // face removed
                    }
                    else if (newFacei < -1)
                    {
                        // face merged
                        elems.append(-newFacei - 2);
                    }
                    else
                    {
                        // valid new face
                        elems.append(newFacei);
                    }
                }
                else
                {
                    // Keep -1 elements
                    elems.append(-1);
                }
            }

            faceList_.transfer(elems);
        }
    }
}


void Foam::probes::movePoints(const polyMesh& mesh)
{
    DebugInformation<< "probes: movePoints" << endl;

    if (fixedLocations_ && &mesh == &mesh_)
    {
        meshUpdated_ = true;
    }
}


// ************************************************************************* //
