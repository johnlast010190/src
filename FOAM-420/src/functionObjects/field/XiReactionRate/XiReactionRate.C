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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "XiReactionRate/XiReactionRate.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(XiReactionRate, 0);
    addToRunTimeSelectionTable(functionObject, XiReactionRate, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::XiReactionRate::XiReactionRate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::XiReactionRate::~XiReactionRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::XiReactionRate::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


bool Foam::functionObjects::XiReactionRate::execute()
{
    return true;
}


bool Foam::functionObjects::XiReactionRate::write()
{
    const volScalarField& b =
        mesh_.lookupObject<volScalarField>("b");

    const volScalarField& Su =
        mesh_.lookupObject<volScalarField>("Su");

    const volScalarField& Xi =
        mesh_.lookupObject<volScalarField>("Xi");

    volScalarField St
    (
        IOobject
        (
            "St",
            time_.timeName(),
            mesh_
        ),
        Xi*Su
    );

    Log << "    Writing turbulent flame-speed field " << St.name()
        << " to " << time_.timeName() << endl;

    St.write();

    volScalarField wdot
    (
        IOobject
        (
            "wdot",
            time_.timeName(),
            mesh_
        ),
        St*Foam::mag(fvc::grad(b))
    );

    Log << "    Writing reaction-rate field " << wdot.name()
        << " to " << time_.timeName() << endl;

    wdot.write();

    return true;
}


// ************************************************************************* //
