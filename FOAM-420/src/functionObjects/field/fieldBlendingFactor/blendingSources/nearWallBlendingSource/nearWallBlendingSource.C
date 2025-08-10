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
#include "fieldBlendingFactor/blendingSources/nearWallBlendingSource/nearWallBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcAverage.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(nearWallBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, nearWallBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nearWallBlendingSource::nearWallBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    nSmooth_(dict.lookupOrDefault<label>("nSmoothIter", 0)),
    fZone_()
{
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
        dimensionedScalar("zi", dimless, 1.),
        zeroGradientFvPatchScalarField::typeName
    );


    forAll(mesh.boundary(), patchI)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchI];

        if (isA<wallPolyPatch>(patch))
        {
            const labelList& meshPoints = patch.meshPoints();
            forAll(meshPoints, mpI)
            {
                label meshPointI = meshPoints[mpI];
                const labelList& pCells = mesh.pointCells()[meshPointI];
                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    cellIndicator[cellI] = 0;
                }
            }
        }
    }

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
                mesh,
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
                    mesh,
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

tmp<surfaceScalarField> nearWallBlendingSource::sourceField()
{
    tmp<surfaceScalarField> zoneIndicator(fZone_);


    return zoneIndicator;
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
