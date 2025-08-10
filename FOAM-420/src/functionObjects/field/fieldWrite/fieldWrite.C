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
    (c) 2010-2012, 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "fieldWrite/fieldWrite.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldWrite, 0);
    addToRunTimeSelectionTable(functionObject, fieldWrite, dictionary);
}
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum<Foam::functionObjects::fieldWrite::spanControls, 5>::names[] =
{
    "startAndEnd",
    "startAndNwrites",
    "start",
    "end",
    "none"
};

const Foam::NamedEnum<Foam::functionObjects::fieldWrite::spanControls, 5>
Foam::functionObjects::fieldWrite::spanControlNames_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldWrite::fieldWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    spanControl_(scNone),
    maxWrites_(label(GREAT)),
    startWrite_(0),
    endWrite_(GREAT),
    fieldNames_(0),
    nWrites_(0),
    outputTime_(false),
    currentValue_(0),
    outputIndex_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldWrite::~fieldWrite()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldWrite::execute()
{
    // Do nothing - only valid on write
    return true;
}

void Foam::functionObjects::fieldWrite::calculate()
{
    outputTime_ = false;

    bool inSpan = false;
    switch(spanControl_)
    {
        case scStartAndEnd:
            if (currentValue_ >= startWrite_ && currentValue_ <= endWrite_)
            {
                inSpan = true;
            }
        break;

        case scStartAndNwrites:
            if (currentValue_ >= startWrite_ && nWrites_ < maxWrites_)
            {
                inSpan = true;
            }
        break;

        case scStart:
            if (currentValue_ >= startWrite_)
            {
                inSpan = true;
            }
        break;

        case scEnd:
            if (currentValue_ <= endWrite_)
            {
                inSpan = true;
            }
        break;

        case scNone:
            inSpan = true;
        break;
    }

    if (inSpan)
    {
        outputTime_ = true;
    }
}

bool Foam::functionObjects::fieldWrite::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    //calculate
    calculate();

    //write
    if (outputTime_)
    {
        nWrites_++;

        if (mesh_.moving())
        {
            lookupObject<pointIOField>(word("points")).write();
        }

        forAll(fieldNames_, fnI)
        {
            if (foundObject<volScalarField>(fieldNames_[fnI]))
            {
                lookupObject<volScalarField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<volVectorField>(fieldNames_[fnI]))
            {
                lookupObject<volVectorField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<volTensorField>(fieldNames_[fnI]))
            {
                lookupObject<volTensorField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<volSymmTensorField>(fieldNames_[fnI]))
            {
                lookupObject<volSymmTensorField>(fieldNames_[fnI]).write();
            }
            else if
            (
                foundObject<volSphericalTensorField>(fieldNames_[fnI])
            )
            {
                lookupObject<volSphericalTensorField>
                    (fieldNames_[fnI]).write();
            }
            else if (foundObject<surfaceScalarField>(fieldNames_[fnI]))
            {
                lookupObject<surfaceScalarField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<surfaceVectorField>(fieldNames_[fnI]))
            {
                lookupObject<surfaceVectorField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<surfaceTensorField>(fieldNames_[fnI]))
            {
                lookupObject<surfaceTensorField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<surfaceSymmTensorField>(fieldNames_[fnI]))
            {
                lookupObject<surfaceSymmTensorField>
                    (fieldNames_[fnI]).write();
            }
            else if
            (
                foundObject<surfaceSphericalTensorField>(fieldNames_[fnI])
            )
            {
                lookupObject<surfaceSphericalTensorField>
                    (fieldNames_[fnI]).write();
            }
            /*
            else if (foundObject<pointScalarField>(fieldNames_[fnI]))
            {
               (lookupObject<pointScalarField>
                    (fieldNames_[fnI])).write();
            }
            else if (foundObject<pointVectorField>(fieldNames_[fnI]))
            {
                lookupObject<pointVectorField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<pointTensorField>(fieldNames_[fnI]))
            {
                lookupObject<pointTensorField>(fieldNames_[fnI]).write();
            }
            else if (foundObject<pointSymmTensorField>(fieldNames_[fnI]))
            {
                lookupObject<pointSymmTensorField>
                    (fieldNames_[fnI]).write();
            }
            else if
            (
                foundObject<pointSphericalTensorField>(fieldNames_[fnI])
            )
            {
                lookupObject<pointSphericalTensorField>
                    (fieldNames_[fnI]).write();
            }*/
            else
            {
                Warning << "No field object named " << fieldNames_[fnI]
                        << " found in object registry." << endl;
            }
        }
    }

    return true;
}


bool Foam::functionObjects::fieldWrite::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    //required entries
    {
        wordList tfns = dict.lookup("fields");
        fieldNames_.setSize(tfns.size());
        fieldNames_ = tfns;
    }


    if (dict.found("spanControl"))
    {
        spanControl_= spanControlNames_.read
        (
            dict.lookup("spanControl")
        );
    }
    else
    {
        spanControl_ = scNone;
    }

    //other entries depend on spanControl_
    switch(spanControl_)
    {
        case scStartAndEnd:
            startWrite_ = readScalar(dict.lookup("start"));
            endWrite_ = readScalar(dict.lookup("end"));
            break;

        case scStartAndNwrites:
            startWrite_ = readScalar(dict.lookup("start"));
            maxWrites_ = readLabel(dict.lookup("maxWrites"));
            break;

        case scStart:
            startWrite_ = readScalar(dict.lookup("start"));
            break;

        case scEnd:
            endWrite_ = readScalar(dict.lookup("end"));
            break;

        case scNone:
            break;

        default:
            FatalError << "Invalid selection of spanControl in fieldWrite."
                       << "Valid options are: StartAndEnd, StartAndNwrites, "
                       << "Start, End and None." << exit(FatalError);
            break;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
