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
    (c) 2011 OpenFOAM Foundation
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionDiffusivity/inverseDistAndVol/inverseDistAndVolDiffusivity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "patchDist/patchDistWave/patchDistWave.H"
#include "patchDist/wallPoint/wallPoint.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseDistAndVolDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseDistAndVolDiffusivity,
        Istream
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseDistAndVolDiffusivity::inverseDistAndVolDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    uniformDiffusivity(mesh, mdData),
    patchNames_(mdData),
    powDist_(readScalar(mdData)),
    powVol_(readScalar(mdData))
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseDistAndVolDiffusivity::~inverseDistAndVolDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::inverseDistAndVolDiffusivity::y() const
{
    const labelHashSet patchSet(mesh().boundaryMesh().patchSet(patchNames_));

    if (patchSet.size())
    {
        tmp<scalarField> tY(new scalarField(mesh().nCells()));
        patchDistWave::wave<wallPoint>(mesh(), patchSet, tY.ref(), false);

        return tY;
    }
    else
    {
        return tmp<scalarField>(new scalarField(mesh().nCells(), 1.0));
    }
}


void Foam::inverseDistAndVolDiffusivity::correct()
{
    volScalarField y_
    (
        IOobject
        (
            "y",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );
    y_.primitiveFieldRef() = y();
    y_.correctBoundaryConditions();

    volScalarField V
    (
        IOobject
        (
            "V",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );

    V.primitiveFieldRef() = mesh().V();
    V.correctBoundaryConditions();

    if (powDist_ != 1. || powVol_ != 1.)
    {
        faceDiffusivity_ = 1.0
            /(pow(fvc::interpolate(V), powVol_) *
            pow(fvc::interpolate(y_), powDist_));
    }
    else
    {
        faceDiffusivity_ = 1.0
            /(fvc::interpolate(V)*fvc::interpolate(y_));
    }
}


// ************************************************************************* //
