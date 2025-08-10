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
    (c) 2011 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingSources/courantNumberBlendingSource/courantNumberBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(courantNumberBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, courantNumberBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

courantNumberBlendingSource::courantNumberBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> courantNumberBlendingSource::sourceField()
{
    tmp<surfaceScalarField> CoNum;

    const surfaceScalarField& phi =
        obr_.lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        CoNum = tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "Cof",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                (
                    mesh_.surfaceInterpolation::deltaCoeffs()
                  * mag(phi)/mesh_.magSf() * mesh_.time().deltaT()
                )
            )
        );
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>(rhoName_);

        CoNum = tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "Cof",
                    mesh_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                (
                    mesh_.surfaceInterpolation::deltaCoeffs()
                  * mag(phi)/(fvc::interpolate(rho)*mesh_.magSf())
                  * mesh_.time().deltaT()
                )
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are incorrect" << nl
            << nl << exit(FatalError);
    }

    return CoNum;
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
