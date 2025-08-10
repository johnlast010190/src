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
    (c) 2018 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "centerOfGravity/centerOfGravity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(centerOfGravity, 0);
    addToRunTimeSelectionTable(functionObject, centerOfGravity, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::centerOfGravity::writeFileHeader(Ostream& os) const
{
    writeCommented(os, "Time");
    writeDelimited(os,"Center of Gravity");
    os << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::centerOfGravity::centerOfGravity
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::centerOfGravity::~centerOfGravity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::centerOfGravity::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    vector CofG = vector::zero;
    vector nsum = vector::zero;
    scalar dsum = scalar(0);

    const scalarField& vol = mesh_.V();
    const vectorField& cc = mesh_.C();

    if (foundObject<volScalarField>(rhoName_))
    {
        // compressible
        const volScalarField& rho =
            lookupObject<volScalarField>(rhoName_);

        forAll(cc, celli)
        {
            scalar cellrho = rho.primitiveField()[celli];
            nsum += cellrho*vol[celli]*cc[celli];
            dsum += cellrho*vol[celli];
        }
    }
    else
    {
        forAll(cc, celli)
        {
            nsum += vol[celli]*cc[celli];
            dsum += vol[celli];
        }
    }

    reduce(nsum, sumOp<vector>());
    reduce(dsum, sumOp<scalar>());

    CofG = nsum/dsum;

    // Write to screen
    Log << "Location of Center Of Gravity: " << CofG <<endl;
    Log << endl;

    writeTime(file());
    file() << CofG;
    file() << endl;

    this->setResult("CofG", CofG);

    return true;
}

bool Foam::functionObjects::centerOfGravity::write()
{
    return true;
}


bool Foam::functionObjects::centerOfGravity::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    return true;
}


// ************************************************************************* //
