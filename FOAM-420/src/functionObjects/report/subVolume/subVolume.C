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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "subVolume/subVolume.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(subVolume, 0);
    addToRunTimeSelectionTable(functionObject, subVolume, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::subVolume::subVolume
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    binary_(dict.lookupOrDefault<bool>("binary", false)),
    writeProcAddressingData_
    (
        dict.lookupOrDefault<bool>("writeProcAddressingData", false)
    ),
    fieldNames_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::subVolume::~subVolume()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::subVolume::execute()
{
    // Do nothing - only valid on write
    return true;
}

void Foam::functionObjects::subVolume::setupSubCase(const word& subName)
{
    //create sub volume dir if it does not exist
    fileName subVolDir;
    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        subVolDir = time_.path()/".."/subName;
    }
    else
    {
        subVolDir = time_.path()/subName;
    }

    // make sub volume folder
    mkDir(subVolDir);
    if (Pstream::parRun())
    {
        Foam::cp
        (
            time_.rootPath()/time_.caseName()/".."/"system",
            subVolDir
        );

        //- create an empty constant folder for paraview readers
        mkDir(subVolDir/word("constant"));
    }
    else
    {
        Foam::cp
        (
            time_.rootPath()/time_.caseName()/"system",
            subVolDir
        );

        //- create an empty constant folder for paraview readers
        mkDir(subVolDir/word("constant"));
    }

    // create sub volume processor directories
    if (Pstream::parRun() && Pstream::master())
    {
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            fileName procDir
            (
                subVolDir/(word("processor") + Foam::name(procI))
            );
            mkDir(procDir);

            //- create an empty constant folder in each proc folder
            mkDir(procDir/word("constant"));
        }
    }

    // Create sub-volume registry
    if (Pstream::parRun())
    {
        subMeshTime_.reset
        (
            new Time
            (
                time_.rootPath()/time_.caseName()/"..",
                subName/fileName
                (word("processor") + Foam::name(Pstream::myProcNo()))
            )
        );
    }
    else
    {
        subMeshTime_.reset
        (
            new Time
            (
                time_.rootPath()/time_.caseName(),
                subName
            )
        );
    }

    // remove function objects from sub case
    {
        IOdictionary controlDict
        (
            IOobject
            (
                "controlDict",
                subMeshTime_().caseSystem(),
                subMeshTime_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if (controlDict.found("functions"))
        {
            controlDict.remove("functions");
            if (binary_)
            {
                controlDict.Foam::regIOobject::writeObject
                (
                    IOstream::BINARY,
                    IOstream::currentVersion,
                    IOstream::COMPRESSED,
                    true
                );
            }
            else
            {
                controlDict.Foam::regIOobject::writeObject
                (
                    IOstream::ASCII,
                    IOstream::currentVersion,
                    IOstream::UNCOMPRESSED,
                    true
                );
            }
        }
    }
}


bool Foam::functionObjects::subVolume::write()
{
    subMeshTime_().setTime(time_.value(), 0);
    if (binary_)
    {
        subMeshTime_().setIOStream
        (
            IOstream::BINARY,
            IOstream::currentVersion,
            IOstream::COMPRESSED
        );
    }
    writeSubset<scalar>(scalarFields_);
    writeSubset<vector>(vectorFields_);
    writeSubset<sphericalTensor>(sphericalTensorFields_);
    writeSubset<symmTensor>(symmTensorFields_);
    writeSubset<tensor>(tensorFields_);

    return true;
}


bool Foam::functionObjects::subVolume::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ":" << nl;

    if (dict.found("fieldNames"))
    {
        WarningInFunction
            << "Old-style of field names  definition (via fieldNames). "
            << " This will no longer be supported in future versions. "
            <<nl
            <<"    Consider using the fields entry. "
            <<endl;
        fieldNames_ = wordList(dict.lookup("fieldNames"));
    }
    if (dict.found("fields"))
    {
        fieldNames_ = wordList(dict.lookup("fields"));
    }
    // Establish type of fields
    findFields<volScalarField>(scalarFields_);
    findFields<volVectorField>(vectorFields_);
    findFields<volSphericalTensorField>(sphericalTensorFields_);
    findFields<volSymmTensorField>(symmTensorFields_);
    findFields<volTensorField>(tensorFields_);

    if (dict.found("writeInterval"))
    {
        writeInterval_ = readLabel(dict.lookup("writeInterval"));
    }
    else
    {
        writeInterval_ = 1;
    }

    setupSubCase(name());

    PtrList<entry> regions(dict.lookup("sets"));

    //- Add all regions to create single sub-setted mesh

    cellSet selectedCellSet
    (
        mesh_,
        "cellSet",
        mesh_.nCells()/10+1  // Reasonable size estimate.
    );

    forAll(regions, regionI)
    {
        const entry& region = regions[regionI];

        autoPtr<topoSetSource> cellSelector =
            topoSetSource::New(region.keyword(), mesh_, region.dict());

        cellSelector->applyToSet
        (
            topoSetSource::ADD,
            selectedCellSet
        );
    }

    label patchI = -1;

    word defaultName("oldInternalFaces");
    if (dict.found("patch"))
    {
        word patchName(dict.lookup("patch"));

        patchI = mesh_.boundaryMesh().findPatchID(patchName);

        if (patchI == -1)
        {
            defaultName = patchName;
        }
        Info<< "Adding exposed internal faces to patch "
            << patchName << endl
            << endl;
    }
    else
    {
        Info<< "Adding exposed internal faces to a patch called"
            << " \"oldInternalFaces\" (created if nessecary)" << endl
            << endl;
    }

    // Create mesh subsetting engine
    subsetMesh_.reset(new fvMeshSubset(mesh_));

    subsetMesh_().setLargeCellSubset
    (
        selectedCellSet,
        patchI,
        true,
        defaultName,
        &subMeshTime_()
    );

    if (binary_)
    {
        subMeshTime_().setIOStream
        (
            IOstream::BINARY,
            IOstream::currentVersion,
            IOstream::COMPRESSED
        );
    }

    if (writeProcAddressingData_)
    {
        subsetMesh_().write();
    }
    else
    {
        subsetMesh_().subMesh().write();
    }

    return true;
}

// ************************************************************************* //
