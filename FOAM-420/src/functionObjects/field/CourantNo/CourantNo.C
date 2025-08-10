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

#include "CourantNo/CourantNo.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(CourantNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        CourantNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::functionObjects::CourantNo::byRho
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


bool Foam::functionObjects::CourantNo::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        const surfaceScalarField& phi =
            lookupObject<surfaceScalarField>(fieldName_);

        tmp<volScalarField::Internal> Coi
        (
            byRho
            (
                (0.5*mesh_.time().deltaT())
               *fvc::surfaceSum(Foam::mag(phi))()()
               /mesh_.V()
            )
        );

        if (foundObject<volScalarField>(resultName_, false))
        {
            volScalarField& Co
            (
                const_cast<volScalarField&>
                (
                    lookupObject<volScalarField>(resultName_)
                )
            );

            Co.ref() = Coi();
            Co.correctBoundaryConditions();
        }
        else
        {
            tmp<volScalarField> tCo
            (
                new volScalarField
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
                    dimensionedScalar("0", dimless, 0.0),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            tCo.ref().ref() = Coi();
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

Foam::functionObjects::CourantNo::CourantNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "phi"),
    rhoName_("rho")
{
    setResultName("Co", "phi");
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::CourantNo::~CourantNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::CourantNo::read(const dictionary& dict)
{
    fieldExpression::read(dict);

    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    return true;
}


// ************************************************************************* //
