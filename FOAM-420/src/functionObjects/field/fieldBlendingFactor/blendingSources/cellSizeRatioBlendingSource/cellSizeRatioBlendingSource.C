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
#include "fieldBlendingFactor/blendingSources/cellSizeRatioBlendingSource/cellSizeRatioBlendingSource.H"
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

defineTypeNameAndDebug(cellSizeRatioBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, cellSizeRatioBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void cellSizeRatioBlendingSource::updateSourceField()
{
    //calculate the ratio of the largest cell volume to the smaller cell
    //volume for each face

    volScalarField cellVol
    (
        IOobject
        (
            "cellVolume",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("vol", dimVolume, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    cellVol.primitiveFieldRef() = mesh_.V();
    cellVol.correctBoundaryConditions();

    string localMax("localMax");
    surfaceScalarField faceMax( fvc::interpolate(cellVol, IStringStream(localMax)()) );

    string localMin("localMin");
    surfaceScalarField faceMin( fvc::interpolate(cellVol, IStringStream(localMin)()) );

    surfaceScalarField::Boundary& faceminbf = faceMin.boundaryFieldRef();
    surfaceScalarField::Boundary& facemaxbf = faceMax.boundaryFieldRef();

    forAll(faceMax.boundaryField(), patchI)
    {
        if (!faceMax.boundaryField()[patchI].coupled())
        {
            facemaxbf[patchI]
                = cellVol.boundaryField()[patchI].patchInternalField();
            faceminbf[patchI]
                = cellVol.boundaryField()[patchI].patchInternalField();
        }
    }

    fZone_.reset((faceMin/faceMax).ptr());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cellSizeRatioBlendingSource::cellSizeRatioBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    fZone_()
{
    updateSourceField();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> cellSizeRatioBlendingSource::sourceField()
{
    if (mesh_.changing())
    {
        updateSourceField();
    }

    tmp<surfaceScalarField> zoneIndicator(fZone_());

    return zoneIndicator;
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
