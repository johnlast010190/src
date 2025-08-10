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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/wallDist/patchDistMethods/Poisson/PoissonPatchDistMethod.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(Poisson, 0);
    addToRunTimeSelectionTable(patchDistMethod, Poisson, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethods::Poisson::Poisson
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs),
    nCorrectors_(dict.lookupOrDefault<label>("nCorrectors", 5)),
    residualTolerance_(dict.lookupOrDefault<scalar>("residualTolerance", 1e-3))
{}


Foam::patchDistMethods::Poisson::Poisson
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    label nCorrectors,
    scalar residualTolerance
)
:
    patchDistMethod(mesh, patchIDs),
    nCorrectors_(nCorrectors),
    residualTolerance_(residualTolerance)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchDistMethods::Poisson::correct(volScalarField& y)
{
    return correct(y, const_cast<volVectorField&>(volVectorField::null()));
}


bool Foam::patchDistMethods::Poisson::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    if (!tyPsi_.valid())
    {
        tyPsi_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "yPsi",
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("yPsi", sqr(dimLength), 0.0),
                y.boundaryFieldRef().types()
            )
        );
    }

    volScalarField& yPsi = tyPsi_.ref();

    //loop to deal with non-orthogonal correction
    for (label nCorr = 0; nCorr < nCorrectors_; nCorr++)
    {
        fvScalarMatrix yPsiEqn
        (
            fvm::laplacian(yPsi) == dimensionedScalar("1", dimless, -1.0)
        );

        scalar initRes = yPsiEqn.solve().initialResidual();
        if (initRes < residualTolerance_)
        {
            break;
        }
    }

    volVectorField gradyPsi(fvc::grad(yPsi));
    volScalarField magGradyPsi("maggradyPsi", mag(gradyPsi));

    y = sqrt(magSqr(gradyPsi) + 2*yPsi) - magGradyPsi;

    // Need to stabilise the y for overset meshes since the holed cells
    // keep the initial value (0.0) so the gradient of that will be
    // zero as well. Turbulence models do not like zero wall distance.
    y.max(SMALL);

    // For overset: enforce smooth y field (yPsi is smooth, magGradyPsi is
    // not)
    mesh_.interpolate(y);


    if (debug && mesh_.time().outputTime()) {y.write();}

    // Cache yPsi if the mesh is moving otherwise delete
    if (!mesh_.changing())
    {
        tyPsi_.clear();
    }

    // Only calculate n if the field is defined
    if (notNull(n))
    {
        n =
           -gradyPsi
           /max
            (
                magGradyPsi,
                dimensionedScalar("smallMagGradyPsi", dimLength, SMALL)
            );

        // For overset: enforce smooth field
        mesh_.interpolate(n);
    }

    return true;
}


// ************************************************************************* //
