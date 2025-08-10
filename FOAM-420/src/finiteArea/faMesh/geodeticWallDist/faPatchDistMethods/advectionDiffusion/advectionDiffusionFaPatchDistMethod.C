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
    (c) 2016 Esi Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "faMesh/geodeticWallDist/faPatchDistMethods/advectionDiffusion/advectionDiffusionFaPatchDistMethod.H"
#include "interpolation/edgeInterpolation/edgeInterpolate.H"
#include "finiteArea/fac/facGrad.H"
#include "finiteArea/fac/facDiv.H"
#include "finiteArea/fam/famDiv.H"
#include "finiteArea/fam/famLaplacian.H"
#include "finiteArea/fam/famSup.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

#include "fields/faPatchFields/basic/fixedValue/fixedValueFaPatchFields.H"
#include "fields/faPatchFields/basic/zeroGradient/zeroGradientFaPatchFields.H"
#include "faMatrices/faMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace faPatchDistMethods
{
    defineTypeNameAndDebug(advectionDiffusion, 0);
    addToRunTimeSelectionTable(faPatchDistMethod, advectionDiffusion, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faPatchDistMethods::advectionDiffusion::advectionDiffusion
(
    const dictionary& dict,
    const faMesh& mesh,
    const labelHashSet& patchIDs
)
:
    faPatchDistMethod(mesh, patchIDs),
    coeffs_(dict.subDict(type() + "Coeffs")),
    pdmPredictor_
    (
        faPatchDistMethod::New
        (
            coeffs_,
            mesh,
            patchIDs
        )
    ),
    epsilon_(coeffs_.lookupOrDefault<scalar>("epsilon", 0.1)),
    tolerance_(coeffs_.lookupOrDefault<scalar>("tolerance", 1e-3)),
    maxIter_(coeffs_.lookupOrDefault<int>("maxIter", 10)),
    predicted_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faPatchDistMethods::advectionDiffusion::correct(areaScalarField& y)
{
    return correct(y, const_cast<areaVectorField&>(areaVectorField::null()));
}


bool Foam::faPatchDistMethods::advectionDiffusion::correct
(
    areaScalarField& y,
    areaVectorField& n
)
{
    if (!predicted_)
    {
        pdmPredictor_->correct(y);
        predicted_ = true;
    }

    areaVectorField ny
    (
        IOobject
        (
            "ny",
            mesh_.time().timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("ny", dimless, vector::zero),
        patchTypes<vector>(mesh_, patchIDs_)
    );

    const faPatchList& patches = mesh_.boundary();
    areaVectorField::Boundary& nybf = ny.boundaryFieldRef();

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();
        nybf[patchi].forceAssign(-patches[patchi].edgeNormals());
    }

    int iter = 0;
    scalar initialResidual = 0;

    do
    {
        ny = fac::grad(y);
        ny /= (mag(ny) + SMALL);

        edgeVectorField nf(fac::interpolate(ny));
        nf /= (mag(nf) + SMALL);

        edgeScalarField yPhi("yPhi", mesh_.Le() & nf);

        faScalarMatrix yEqn
        (
            fam::div(yPhi, y)
          - fam::Sp(fac::div(yPhi), y)
          - epsilon_*y*fam::laplacian(y)
         ==
            dimensionedScalar("1", dimless, 1.0)
        );

        yEqn.relax();
        initialResidual = yEqn.solve().initialResidual();

    } while (initialResidual > tolerance_ && ++iter < maxIter_);

    // Only calculate n if the field is defined
    if (notNull(n))
    {
        n = -ny;
    }

    return true;
}


// ************************************************************************* //
