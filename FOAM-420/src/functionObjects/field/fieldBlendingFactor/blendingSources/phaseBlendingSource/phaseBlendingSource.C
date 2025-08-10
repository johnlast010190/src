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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingSources/phaseBlendingSource/phaseBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "incompressibleTwoPhaseMixture/incompressibleTwoPhaseMixture.H"
#include "interpolation/surfaceInterpolation/schemes/localMax/localMax.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(phaseBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, phaseBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

phaseBlendingSource::phaseBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict)
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> phaseBlendingSource::sourceField()
{

    tmp<surfaceScalarField> Interface;

    if
    (
        mesh_.foundObject<incompressibleTwoPhaseMixture>
        ("transportProperties")
    )
    {
        const incompressibleTwoPhaseMixture& mixture =
            obr_.lookupObject<incompressibleTwoPhaseMixture>
            ("transportProperties");

        scalar df = 0.1;
        scalar dfi = 1 - df;

        dimensionedScalar rhoMin
        (
            "rhoMin",
            mixture.rho1().dimensions(),
            min
            (
                mixture.rho1().value(),
                mixture.rho2().value()
            )
        );
        dimensionedScalar rhoMax
        (
            "rhoMax",
            mixture.rho1().dimensions(),
            max
            (
                mixture.rho1().value(),
                mixture.rho2().value()
            )
        );

        const volScalarField& rho
        (
            obr_.lookupObject<volScalarField>("rho")
        );


        volScalarField intc
        (
            "intc",
            pos0(rhoMax*dfi - rho)
            *pos0(rho - rhoMin*(1+df))
        );

        Interface = tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phaseInterface",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                localMax<scalar>(mesh_).interpolate(intc)
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "This blending source is only valid in conjunction with "
            << "interFoam. "
            << nl << exit(FatalError);
    }

    return Interface;
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
