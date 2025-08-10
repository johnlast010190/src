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
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::autoPtr<typename Foam::LduMatrix<Type, DType, LUType>::smoother>
Foam::LduMatrix<Type, DType, LUType>::smoother::New
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& smootherDict
)
{
    word smootherName = smootherDict.lookup("smoother");

    if (matrix.symmetric())
    {
        typename symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTable_().find(smootherName);

        if (constructorIter == symMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction(smootherDict)
                << "Unknown symmetric matrix smoother " << smootherName
                << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            constructorIter->second
            (
                fieldName,
                matrix
            )
        );
    }
    else if (matrix.asymmetric())
    {
        typename asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTable_().find(smootherName);

        if (constructorIter == asymMatrixConstructorTable_().end())
        {
            FatalIOErrorInFunction(smootherDict)
                << "Unknown asymmetric matrix smoother " << smootherName
                << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            constructorIter->second
            (
                fieldName,
                matrix
            )
        );
    }
    else
    {
        FatalIOErrorInFunction(smootherDict)
            << "cannot solve incomplete matrix, no off-diagonal coefficients"
            << exit(FatalIOError);

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            nullptr
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::smoother::smoother
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix
)
:
    fieldName_(fieldName),
    matrix_(matrix)
{}


// ************************************************************************* //
