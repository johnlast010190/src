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
    (c) 2011-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "fieldBlendingFactor/blendingSources/zonalDESBlendingSource/zonalDESBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcAverage.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(zonalDESBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, zonalDESBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

zonalDESBlendingSource::zonalDESBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    CDES_(dict.lookupOrDefault<scalar>("CDES", 0.65)),
    nSmooth_(dict.lookupOrDefault<label>("nSmoothIter", 0)),
    fZone_()
{
    const volScalarField& y(wallDist::New(mesh).y());
    scalarField deltaV(mesh.nCells(), 0.0);

    label nD = mesh.nGeometricD();

    if (nD == 3)
    {
        deltaV = CDES_*pow(mesh.V(), 1.0/3.0);
    }
    else if (nD == 2)
    {
        const Vector<label>& directions = mesh.geometricD();

        scalar thickness = 0.0;
        for (direction dir=0; dir<directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh.bounds().span()[dir];
                break;
            }
        }

        deltaV = CDES_*sqrt(mesh.V()/thickness);
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, not applicable"
            << exit(FatalError);
    }

    // zet to 0 for near-wall cell and 1 otherwise
    volScalarField cellIndicator
    (
        IOobject
        (
            "zoneIndicator",
            mesh.time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zi", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    cellIndicator.primitiveFieldRef() = pos0(y.internalField()-deltaV);
    cellIndicator.correctBoundaryConditions();

    cellIndicator.write();

    fZone_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "zoneFaceIndicator",
                mesh.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::interpolate(cellIndicator)
        )
    );

    for (label i = 0; i < nSmooth_; i++)
    {
        fZone_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "zoneFaceIndicator",
                    mesh.time().timeName(),
                    obr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(cellIndicator)
            )
        );
        cellIndicator.primitiveFieldRef()=
            (fvc::average(fZone_()))->internalField();
        cellIndicator.correctBoundaryConditions();
    }


}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> zonalDESBlendingSource::sourceField()
{
    tmp<surfaceScalarField> zoneIndicator(fZone_);


    return zoneIndicator;
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
