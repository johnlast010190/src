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
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "patchWrite/patchWrite.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(patchWrite, 0);
    addToRunTimeSelectionTable(functionObject, patchWrite, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::patchWrite::patchWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writePatchNames_(0),
    writePatchIDs_(0),
    fieldNames_(0),
    nearCellValue_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::patchWrite::~patchWrite()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::patchWrite::execute()
{
    // Do nothing - only valid on write
    return true;
}


void Foam::functionObjects::patchWrite::calculate()
{
}

bool Foam::functionObjects::patchWrite::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    //calculate

    forAll(writePatchIDs_, wpiI)
    {
        if (writePatchIDs_[wpiI] != -1)
        {
            label patchID = writePatchIDs_[wpiI];

            word postfix = "";

            if (nearCellValue_)
            {
                postfix = "_nw";
            }

            mkDir(obr_.time().path()/writePatchNames_[wpiI]);

            fileName writeDir
            (
                obr_.time().path()/writePatchNames_[wpiI]/obr_.time().timeName()
            );

            mkDir(writeDir);

            forAll(fieldNames_, fnI)
            {
                OFstream patchFile
                (
                    writeDir/fieldNames_[fnI]+postfix,
                    IOstream::ASCII,
                    IOstream::currentVersion,
                    IOstream::COMPRESSED
                );

                if (foundObject<volScalarField>(fieldNames_[fnI]))
                {
                    const volScalarField& cf
                        = lookupObject<volScalarField>(fieldNames_[fnI]);
                    const scalarField& pf = cf.boundaryField()[patchID];
                    patchFile << pf;

                }
                else if (foundObject<volVectorField>(fieldNames_[fnI]))
                {
                    const volScalarField& cf
                        = lookupObject<volScalarField>(fieldNames_[fnI]);
                    const scalarField& pf = cf.boundaryField()[patchID];
                    patchFile << pf;
                }
                else if (foundObject<volTensorField>(fieldNames_[fnI]))
                {
                    const volScalarField& cf
                        = lookupObject<volScalarField>(fieldNames_[fnI]);
                    const scalarField& pf = cf.boundaryField()[patchID];
                    patchFile << pf;
                }
                else
                {
                    WarningInFunction
                        << "patchWrite:" << name()
                        << " cannot find a volField called "
                        << fieldNames_[fnI]
                        << " in the object registry. Skipping." << endl;
                }
            }
        }
    }
    return true;
}


bool Foam::functionObjects::patchWrite::read(const dictionary& dict)
{
    Log << type() << " " << name() <<  " read:" << nl;

    writePatchNames_ = wordList(dict.lookup("patches"));

    writePatchIDs_.setSize(writePatchNames_.size());

    forAll(writePatchNames_, wpnI)
    {
        writePatchIDs_[wpnI]
            = mesh_.boundaryMesh().findPatchID(writePatchNames_[wpnI]);

        if (writePatchIDs_[wpnI] == -1)
        {
            WarningInFunction
                << "Could not find patch with name "
                << writePatchNames_[wpnI] << " in boundary list." << nl
                << "Entry will be ignored by patchWrite:"<< name() << "."
                << endl;
        }
    }

    fieldNames_ = wordList(dict.lookup("fields"));

    nearCellValue_ = false;

    if (dict.found("nearCellValues"))
    {
        nearCellValue_ = readBool(dict.lookup("nearCellValues"));
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
