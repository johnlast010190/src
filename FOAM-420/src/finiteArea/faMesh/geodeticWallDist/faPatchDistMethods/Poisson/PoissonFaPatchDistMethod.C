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

#include "finiteArea/fa/fa.H"
#include "faMesh/geodeticWallDist/faPatchDistMethods/Poisson/PoissonFaPatchDistMethod.H"
#include "finiteArea/fac/facGrad.H"
#include "finiteArea/fam/famLaplacian.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "faMatrices/faMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace faPatchDistMethods
{
    defineTypeNameAndDebug(Poisson, 0);
    addToRunTimeSelectionTable(faPatchDistMethod, Poisson, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faPatchDistMethods::Poisson::Poisson
(
    const dictionary& dict,
    const faMesh& mesh,
    const labelHashSet& patchIDs
)
:
    faPatchDistMethod(mesh, patchIDs)
{}


Foam::faPatchDistMethods::Poisson::Poisson
(
    const faMesh& mesh,
    const labelHashSet& patchIDs
)
:
    faPatchDistMethod(mesh, patchIDs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faPatchDistMethods::Poisson::correct(areaScalarField& y)
{
    return correct(y, const_cast<areaVectorField&>(areaVectorField::null()));
}


bool Foam::faPatchDistMethods::Poisson::correct
(
    areaScalarField& y,
    areaVectorField& n
)
{
    if (!tyPsi_.valid())
    {
        tyPsi_ = tmp<areaScalarField>
        (
            new areaScalarField
            (
                IOobject
                (
                    "yPsi",
                    mesh_.time().timeName(),
                    mesh_.thisDb()
                ),
                mesh_,
                dimensionedScalar("yPsi", sqr(dimLength), 0.0),
                y.boundaryField().types()
            )
        );
    }

    areaScalarField& yPsi = tyPsi_.ref();

       solve(fam::laplacian(yPsi) == dimensionedScalar("1", dimless, -1.0));

    areaVectorField gradyPsi(fac::grad(yPsi));
    areaScalarField magGradyPsi(mag(gradyPsi));

    y = sqrt(magSqr(gradyPsi) + 2*yPsi) - magGradyPsi;

    // Cache yPsi if the mesh is moving otherwise delete
//    if (!mesh_.changing())
//    {
        tyPsi_.clear();
//    }

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
    }

    return true;
}


// ************************************************************************* //
