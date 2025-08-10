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
    (c) 2004-6 H. Jasak All rights reserved

Description
    Solver for diagonal matrices.

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduSolvers/BlockDiagonal/BlockDiagonalSolver.H"
#include "matrices/blockLduMatrix/BlockLduSolvers/BlockLduSolver/BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
template<class Type>
Foam::BlockDiagonalSolver<Type>::BlockDiagonalSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<Type>(fieldName, matrix, dict)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockSolverPerformance<Type>
Foam::BlockDiagonalSolver<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b
)
{
    FatalErrorInFunction
        << "Not supported"
        << abort(FatalError);

    //CoeffField<Type> dD = inv(this->matrix_.diag());

    //multiply(x, b, dD);

    return BlockSolverPerformance<Type>
    (
        this->typeName,
        this->fieldName(),
        pTraits<Type>::zero,
        pTraits<Type>::zero,
        0,
        true,
        false
    );
}


// ************************************************************************* //
