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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "PecletNo/PecletNo.H"
#include "turbulenceModel.H"
// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(PecletNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        PecletNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::functionObjects::PecletNo::rhoScale
(
    const surfaceScalarField& phi
) const
{
    if (phi.dimensions() == dimMass/dimTime)
    {
        return phi/fvc::interpolate(lookupObject<volScalarField>(rhoName_));
    }
    else
    {
        return phi;
    }
}


bool Foam::functionObjects::PecletNo::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        tmp<volScalarField> nuEff;
        if (mesh_.foundObject<turbulenceModel>(turbulenceModel::propertiesName))
        {
            const turbulenceModel& model =
                lookupObject<turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );

            nuEff = model.nuEff();
        }
        else if (mesh_.foundObject<dictionary>("transportProperties"))
        {
            const dictionary& model =
                mesh_.lookupObject<dictionary>("transportProperties");

            nuEff =
                tmp<volScalarField>
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "nuEff",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar("nu", dimViscosity/dimDensity, model)
                    )
                );
        }
        else
        {
            FatalErrorInFunction
                << "Unable to determine the viscosity"
                << exit(FatalError);
        }


        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(fieldName_);

        return store
        (
            resultName_,
            Foam::mag(rhoScale(phi))
           /(
                mesh_.magSf()
               *mesh_.surfaceInterpolation::deltaCoeffs()
               *fvc::interpolate(nuEff)
            )
        );
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNo::PecletNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "phi"),
    rhoName_("rho")
{
    setResultName("Pe", "phi");
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNo::~PecletNo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::PecletNo::read(const dictionary& dict)
{
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    return true;
}


// ************************************************************************* //
