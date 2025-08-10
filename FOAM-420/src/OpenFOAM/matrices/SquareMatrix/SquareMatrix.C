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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "matrices/SquareMatrix/SquareMatrix.H"
#include "primitives/ints/lists/labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::detDecomposed
(
    const SquareMatrix<Type>& matrix,
    const label sign
)
{
    Type diagProduct = pTraits<Type>::one;

    for (label i=0; i<matrix.m(); i++)
    {
        diagProduct *= matrix(i, i);
    }

    return sign*diagProduct;
}


template<class Type>
Foam::scalar Foam::det(const SquareMatrix<Type>& matrix)
{
    SquareMatrix<Type> matrixTmp = matrix;

    labelList pivotIndices(matrix.m());
    label sign;
    LUDecompose(matrixTmp, pivotIndices, sign);

    return detDecomposed(matrixTmp, sign);
}


template<class Type>
Foam::scalar Foam::det(SquareMatrix<Type>& matrix)
{
    labelList pivotIndices(matrix.m());
    label sign;
    LUDecompose(matrix, pivotIndices, sign);

    return detDecomposed(matrix, sign);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
