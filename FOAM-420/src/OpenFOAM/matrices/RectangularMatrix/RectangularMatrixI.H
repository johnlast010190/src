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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix()
:
    Matrix<RectangularMatrix<Type>, Type>()
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const label m,
    const label n
)
:
    Matrix<RectangularMatrix<Type>, Type>(m, n)
{}


template<class Type>
template<class MatrixType>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const ConstMatrixBlock<MatrixType>& block
)
:
    Matrix<RectangularMatrix<Type>, Type>(block)
{}


template<class Type>
template<class MatrixType>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const MatrixBlock<MatrixType>& block
)
:
    Matrix<RectangularMatrix<Type>, Type>(block)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const label m,
    const label n,
    const zero
)
:
    Matrix<RectangularMatrix<Type>, Type>(m, n, Zero)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const label m,
    const label n,
    const Type& t
)
:
    Matrix<RectangularMatrix<Type>, Type>(m, n, t)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix
(
    const SquareMatrix<Type>& SM
)
:
    Matrix<RectangularMatrix<Type>, Type>(SM)
{}


template<class Type>
inline Foam::RectangularMatrix<Type>::RectangularMatrix(Istream& is)
:
    Matrix<RectangularMatrix<Type>, Type>(is)
{}


template<class Type>
inline Foam::autoPtr<Foam::RectangularMatrix<Type>>
Foam::RectangularMatrix<Type>::clone() const
{
    return autoPtr<RectangularMatrix<Type>>
    (
        new RectangularMatrix<Type>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::RectangularMatrix<Type>::operator=(const zero)
{
    Matrix<RectangularMatrix<Type>, Type>::operator=(Zero);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
inline Foam::RectangularMatrix<Type> outer
(
    const Field<Type>& f1,
    const Field<Type>& f2
)
{
    RectangularMatrix<Type> f1f2T(f1.size(), f2.size());

    for (label i=0; i<f1f2T.m(); i++)
    {
        for (label j=0; j<f1f2T.n(); j++)
        {
            f1f2T(i, j) = f1[i]*f2[j];
        }
    }

    return f1f2T;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
