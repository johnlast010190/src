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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/lduMatrix/solvers/GAMG/GAMGSolver.H"
#include "primitives/Vector2D/vector2D/vector2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::scale
(
    scalarField& field,
    scalarField& Acf,
    const lduMatrix& A,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const scalarField& source,
    const direction cmpt
) const
{
    A.Amul
    (
        Acf,
        field,
        interfaceLevelBouCoeffs,
        interfaceLevel,
        cmpt
    );


    const label nCells = field.size();
    scalar* __restrict__ fieldPtr = field.begin();
    const scalar* const __restrict__ sourcePtr = source.begin();
    const scalar* const __restrict__ AcfPtr = Acf.begin();


    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    for (label i=0; i<nCells; i++)
    {
        scalingFactorNum += sourcePtr[i]*fieldPtr[i];
        scalingFactorDenom += AcfPtr[i]*fieldPtr[i];
    }

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    A.mesh().reduce(scalingVector, sumOp<vector2D>());

    const scalar sf = scalingVector.x()/stabilise(scalingVector.y(), VSMALL);

    if (debug >= 2)
    {
        Pout<< sf << " ";
    }

    const scalarField& D = A.diag();
    const scalar* const __restrict__ DPtr = D.begin();

    for (label i=0; i<nCells; i++)
    {
        fieldPtr[i] = sf*fieldPtr[i] + (sourcePtr[i] - sf*AcfPtr[i])/DPtr[i];
    }
}


// ************************************************************************* //
