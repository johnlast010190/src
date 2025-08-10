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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/dictionary/dictionary.H"
#include "residualBlendingSource.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "finiteVolume/fvc/fvcGrad.H"
// Use surfaceInterpolate from finiteVolume, not this library!
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace blendingSources
{

defineTypeNameAndDebug(residualBlendingSource, 0);
addToRunTimeSelectionTable
(blendingSource, residualBlendingSource, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

residualBlendingSource::residualBlendingSource
(
    const objectRegistry& obr,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    blendingSource(obr, mesh, dict),
    residualName_(dict.lookup("residual")),
    error_(mesh_.nCells(), Zero),
    errorIntegral_(mesh_.nCells(), Zero),
    oldError_(mesh_.nCells(), Zero),
    oldErrorIntegral_(mesh_.nCells(), Zero),
    P_(dict.lookupOrDefault<scalar>("P", 3)),
    I_(dict.lookupOrDefault<scalar>("I", 0.0)),
    D_(dict.lookupOrDefault<scalar>("D", 0.25))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<surfaceScalarField> residualBlendingSource::sourceField()
{
    const scalarField& residual = obr_.lookupObject<scalarField>(residualName_);

    scalar meanRes = gAverage(mag(residual)()) + VSMALL;

    oldError_ = error_;
    oldErrorIntegral_ = errorIntegral_;
    error_ = mag(meanRes - mag(residual));
    errorIntegral_ = oldErrorIntegral_ + 0.5*(error_ + oldError_);
    const scalarField errorDifferential(error_ - oldError_);

    const scalarField factorList
    (
        + P_*error_
        + I_*errorIntegral_
        + D_*errorDifferential
    );

    tmp<volScalarField> tres;

    tres = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "residualBlending",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("res", dimless, Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tres.ref().primitiveFieldRef() = factorList;
    tres.ref().correctBoundaryConditions();

    // update blending function
    ramp_->update(meanRes);

    return fvc::interpolate(tres());
}

// ************************************************************************* //

} // End namespace blendingSources
} // End namespace Foam

// ************************************************************************* //
