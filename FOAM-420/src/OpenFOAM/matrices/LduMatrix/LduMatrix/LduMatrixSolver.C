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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"
#include "matrices/LduMatrix/Solvers/DiagonalSolver/DiagonalSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::autoPtr<typename Foam::LduMatrix<Type, DType, LUType>::solver>
Foam::LduMatrix<Type, DType, LUType>::solver::New
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
{
    word solverName = solverDict.lookup("solver");

    if (matrix.diagonal())
    {
        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
        (
            new DiagonalSolver<Type, DType, LUType>
            (
                fieldName,
                matrix,
                solverDict
            )
        );
    }
    else if (matrix.symmetric())
    {
        typename symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTable_().find(solverName);

        if (constructorIter == symMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction(solverDict)
                << "Unknown symmetric matrix solver " << solverName
                << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
        (
            constructorIter->second
            (
                fieldName,
                matrix,
                solverDict
            )
        );
    }
    else if (matrix.asymmetric())
    {
        typename asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTable_().find(solverName);

        if (constructorIter == asymMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction(solverDict)
                << "Unknown asymmetric matrix solver " << solverName
                << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
        (
            constructorIter->second
            (
                fieldName,
                matrix,
                solverDict
            )
        );
    }
    else
    {
        FatalIOErrorInFunction(solverDict)
            << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
        (
            nullptr
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::solver::solver
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
:
    fieldName_(fieldName),
    matrix_(matrix),

    controlDict_(solverDict),

    maxIter_(defaultMaxIter_),
    minIter_(0),
    tolerance_(1e-6*pTraits<Type>::one),
    relTol_(Zero)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::solver::readControls()
{
    readControl(controlDict_, maxIter_, "maxIter");
    readControl(controlDict_, minIter_, "minIter");
    readControl(controlDict_, tolerance_, "tolerance");
    readControl(controlDict_, relTol_, "relTol");
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::solver::read
(
    const dictionary& solverDict
)
{
    controlDict_ = solverDict;
    readControls();
}


template<class Type, class DType, class LUType>
Type Foam::LduMatrix<Type, DType, LUType>::solver::normFactor
(
    const Field<Type>& psi,
    const Field<Type>& Apsi,
    Field<Type>& tmpField
) const
{
    // --- Calculate A dot reference value of psi
    matrix_.sumA(tmpField);
    cmptMultiply(tmpField, tmpField, gAverage(psi));

    return stabilise
    (
        gSum(cmptMag(Apsi - tmpField) + cmptMag(matrix_.source() - tmpField)),
        SolverPerformance<Type>::small_
    );

    // At convergence this simpler method is equivalent to the above
    // return stabilise(2*gSumCmptMag(matrix_.source()), matrix_.small_);
}


// ************************************************************************* //
