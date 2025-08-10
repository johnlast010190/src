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
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "maxFaceDelta.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(maxFaceDelta, 0);
addToRunTimeSelectionTable(LESdelta, maxFaceDelta, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void maxFaceDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();
    label nD = mesh.nGeometricD();

    //zero delta before re-assignment
    delta_.primitiveFieldRef() = 0.0;

    const cellList& cells = mesh.cells();
    const vectorField& Sf(mesh.faceAreas());
    const scalarField& Vol(mesh.V());

    forAll(cells,cellI)
    {
        scalar& deltaMax = delta_[cellI];
        const labelList& cFaces = mesh.cells()[cellI];
        const point& C = mesh.cellCentres()[cellI];

        forAll(cFaces, cFaceI)
        {
            const point& Cf = mesh.faceCentres()[cFaces[cFaceI]];
            vector faceDir = Cf - C;
            faceDir /= mag(faceDir);

            scalar Snorm = 0;

            //second loop to sum projected face areas
            forAll(cFaces, cFaceII)
            {
                Snorm += mag(faceDir & Sf[cFaces[cFaceII]]);
            }
            deltaMax = max(deltaMax, 2*Vol[cellI]/Snorm);
        }
        deltaMax *= deltaCoeff_;
    }

    //update processor boundaries
    delta_.correctBoundaryConditions();

    if (nD == 3)
    {
        //ok
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable\n"
            << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

maxFaceDelta::maxFaceDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dd
)
:
    LESdelta(name, turbulence),
    deltaCoeff_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<scalar>("deltaCoeff", 1.0)
    ),
    recalcDelta_
    (
        dd.subDict(type() + "Coeffs").lookupOrDefault<bool>
        (
            "recomputeDelta",
            true
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void maxFaceDelta::read(const dictionary& dd)
{
    deltaCoeff_
        = dd.subDict(type() + "Coeffs")
          .lookupOrDefault<scalar>("deltaCoeff", 1.0);

    recalcDelta_ = dd.subDict(type() + "Coeffs").lookupOrDefault<bool>
            (
                "recomputeDelta",
                true
            );

    calcDelta();
}


void maxFaceDelta::correct()
{
    if (turbulenceModel_.mesh().changing() && recalcDelta_)
    {
        calcDelta();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
