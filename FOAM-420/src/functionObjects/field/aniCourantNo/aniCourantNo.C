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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "aniCourantNo/aniCourantNo.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(aniCourantNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        aniCourantNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::functionObjects::aniCourantNo::byRho
(
    const tmp<volScalarField::Internal>& Co
) const
{
    if (Co().dimensions() == dimDensity)
    {
        return Co/lookupObject<volScalarField>(rhoName_);
    }
    else
    {
        return Co;
    }
}


bool Foam::functionObjects::aniCourantNo::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        vectorField minC(mesh_.nCells(), vector::zero);
        vectorField maxC(mesh_.nCells(), vector::zero);

        forAll(mesh_.cells(), celli)
        {
            const pointField cp
            (
                mesh_.cells()[celli].points(mesh_.faces(), mesh_.points())
            );
            minC[celli] = min(cp);
            maxC[celli] = max(cp);
        }

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        vectorField Coi(mesh_.nCells(), vector::zero);
        forAll(Coi, celli)
        {
            vector dx = maxC[celli] - minC[celli];
            Coi[celli].x() = Foam::mag(U[celli].x())/dx.x();
            Coi[celli].y() = Foam::mag(U[celli].y())/dx.y();
            Coi[celli].z() = Foam::mag(U[celli].z())/dx.z();
            Coi[celli] *= mesh_.time().deltaT().value();
        }

        if (foundObject<volVectorField>(resultName_, false))
        {
            volVectorField& Co = lookupObjectRef<volVectorField>(resultName_);
            Co.primitiveFieldRef() = Coi;
            Co.correctBoundaryConditions();
        }
        else
        {
            tmp<volVectorField> tCo
            (
                new volVectorField
                (
                    IOobject
                    (
                        resultName_,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("0", dimless, vector(0, 0, 0)),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            tCo->primitiveFieldRef() = Coi;
            tCo.ref().correctBoundaryConditions();
            mesh_.objectRegistry::store(tCo.ptr());
        }
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::aniCourantNo::aniCourantNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "phi"),
    rhoName_("rho"),
    UName_("U")
{
    setResultName("aniCo", "phi");
    aniCourantNo::read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::aniCourantNo::~aniCourantNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::aniCourantNo::read(const dictionary& dict)
{
    fieldExpression::read(dict);

    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    UName_ = dict.lookupOrDefault<word>("UName", "U");

    return true;
}


// ************************************************************************* //
