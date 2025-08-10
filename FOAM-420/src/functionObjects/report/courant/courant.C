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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "courant/courant.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "finiteVolume/fvc/fvc.H"
#include "fieldInstance/fieldInstance.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(courant, 0);
    addToRunTimeSelectionTable(functionObject, courant, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::courant::writeFileHeader(Ostream& os) const
{
    writeCommented(os,"Time MaxCo at location MeanCo");
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::courant::courant
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

Foam::functionObjects::courant::~courant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::courant::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    fieldInstance coLocation(-1);

    autoPtr<surfaceScalarField> CoPtr;
    const auto& phi =
        lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        // compressible
        const auto& rho =
            lookupObject<volScalarField>(rhoName_);

        CoPtr.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "Co",
                    obr_.time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (
                    phi.mesh().surfaceInterpolation::deltaCoeffs()
                    * (mag(phi)/(fvc::interpolate(rho)*phi.mesh().magSf()))
                    * obr_.time().deltaT()
                )
             )
         );
    }
    else if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        // incompressible
        CoPtr.set
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "Co",
                    obr_.time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (
                    phi.mesh().surfaceInterpolation::deltaCoeffs()
                    * (mag(phi)/phi.mesh().magSf())
                    * obr_.time().deltaT()
                )
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);
    }

    const surfaceScalarField& co=CoPtr();
    const surfaceVectorField& cf=co.mesh().Cf();

    forAll(co,faceI)
    {
        if (co[faceI] > coLocation.value())
        {
            coLocation = fieldInstance(co[faceI],cf[faceI]);
        }
    }

    forAll(co.boundaryField(),patchI)
    {
        forAll(co.boundaryField()[patchI], faceI)
        {
            if (co.boundaryField()[patchI][faceI] > coLocation.value())
            {
                coLocation = fieldInstance
                (
                    co.boundaryField()[patchI][faceI],
                    cf.boundaryField()[patchI][faceI]
                );
            }
        }
    }

    // Synchronize max courant number on processors
    reduce(coLocation, maxFIOp());

    scalar meanCo = gAverage(co);

    // Write to screen
    Log << "Maximum Courant: " << coLocation.value() << " at "
        << coLocation.position()
        << ", mean Courant: " << meanCo << endl;
    Log << endl;

    file() << obr_.time().timeName() << " ";
    file() << coLocation.value() << " ";
    file() << coLocation.position() << " ";
    file() << meanCo << endl;

    return true;
}

bool Foam::functionObjects::courant::write()
{
    return true;
}


bool Foam::functionObjects::courant::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    // Setup phi and rho
    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }
    else
    {
        phiName_="phi";
    }

    if (dict.found("rho"))
    {
        dict.lookup("rho") >> rhoName_;
    }
    else
    {
        rhoName_="rho";
    }

    return true;
}


// ************************************************************************* //
