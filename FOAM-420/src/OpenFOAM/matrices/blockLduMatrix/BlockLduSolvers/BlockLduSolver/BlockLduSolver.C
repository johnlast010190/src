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

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduMatrix/blockLduMatrices.H"
#include "matrices/blockLduMatrix/BlockLduSolvers/BlockDiagonal/BlockDiagonalSolver.H"
#include "db/IOstreams/IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::scalar Foam::BlockLduSolver<Type>::great_ = 1.0e+20;

template<class Type>
const Foam::scalar Foam::BlockLduSolver<Type>::small_ = 1.0e-20;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class T>
void Foam::BlockLduSolver<Type>::readControl
(
    const dictionary& dict,
    T& control,
    const word& controlName
)
{
    if (dict.found(controlName))
    {
        dict.lookup(controlName) >> control;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Cannot read control " << controlName
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockLduSolver<Type>::BlockLduSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix
)
:
    fieldName_(fieldName),
    dict_(),
    matrix_(matrix)
{}


template<class Type>
Foam::BlockLduSolver<Type>::BlockLduSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    fieldName_(fieldName),
    dict_(dict),
    matrix_(matrix)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<typename Foam::BlockLduSolver<Type>>
Foam::BlockLduSolver<Type>::New
(
    const word& fieldName,
    BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
{
    word solverName(dict.lookup("solver"));

    if (matrix.diagonal())
    {
        return autoPtr<BlockLduSolver<Type>>
        (
            new BlockDiagonalSolver<Type>(fieldName, matrix, dict)
        );
    }
    else if (matrix.symmetric())
    {
        typename symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTable_().find(solverName);

        if (constructorIter == symMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown symmetric matrix solver " << solverName
                << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<BlockLduSolver<Type>>
        (
            constructorIter->second
            (
                fieldName,
                matrix,
                dict
            )
        );
    }
    else if (matrix.asymmetric())
    {
        if (asymMatrixConstructorTable_().empty())
        {
            FatalErrorInFunction
                << "No Block solvers found. "
                << "Check linked libraries or compilation build"
                << abort(FatalError);
        }

        typename asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTable_().find(solverName);

        if (constructorIter == asymMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown asymmetric matrix solver " << solverName
                << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<BlockLduSolver<Type>>
        (
            constructorIter->second
            (
                fieldName,
                matrix,
                dict
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalError);

        return autoPtr<BlockLduSolver<Type>>(nullptr);
    }
}


// ************************************************************************* //
