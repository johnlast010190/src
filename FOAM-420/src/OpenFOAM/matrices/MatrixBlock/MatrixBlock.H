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
    (c) 2016 OpenFOAM Foundation

Class
    Foam::MatrixBlock

Description
    A templated block of an (m x n) matrix of type \<MatrixType\>.

        Foam::ConstMatrixBlock: block of a const matrix
        Foam::MatrixBlock: block of a non-const matrix

    The block may be assigned to a block of another matrix or to a VectorSpace
    or MatrixSpace e.g. \c tensor.  Conversion of a column block to a \c
    Field<T> is also provide.

SourceFiles
    MatrixBlock.C
    MatrixBlockI.H

\*---------------------------------------------------------------------------*/

#ifndef MatrixBlock_H
#define MatrixBlock_H

#include "matrices/Matrix/Matrix.H"
#include "primitives/MatrixSpace/MatrixSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ConstMatrixBlock Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class ConstMatrixBlock
{
    // Private data

        //- Const reference to the parent matrix
        const MatrixType& matrix_;

        // Block size
        const label mRows_;
        const label nCols_;

        // Block location in parent matrix
        const label rowStart_;
        const label colStart_;

public:

    typedef typename MatrixType::cmptType cmptType;

    // Constructors

        //- Construct block for matrix, size and location
        inline ConstMatrixBlock
        (
            const MatrixType& matrix,
            const label m,
            const label n,
            const label mStart,
            const label nStart
        );


    // Member Functions

        //- Return the number of rows in the block
        inline label m() const;

        //- Return the number of columns in the block
        inline label n() const;

        //- (i, j) const element access operator
        inline const cmptType& operator()
        (
            const label i,
            const label j
        ) const;

        //- Convert a column of a matrix to a Field
        operator Field<cmptType>() const;
};


/*---------------------------------------------------------------------------*\
                          Class MatrixBlock Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class MatrixBlock
{
    // Private data

        //- Reference to the parent matrix
        MatrixType& matrix_;

        // Block size
        const label mRows_;
        const label nCols_;

        // Block location in parent matrix
        const label rowStart_;
        const label colStart_;

public:

    typedef typename MatrixType::cmptType cmptType;

    // Constructors

        //- Construct block for matrix, size and location
        inline MatrixBlock
        (
            MatrixType& matrix,
            const label m,
            const label n,
            const label mStart,
            const label nStart
        );


    // Member Functions

        //- Return the number of rows in the block
        inline label m() const;

        //- Return the number of columns in the block
        inline label n() const;

        //- (i, j) const element access operator
        inline const cmptType& operator()
        (
            const label i,
            const label j
        ) const;

        //- (i, j) element access operator
        inline cmptType& operator()(const label i, const label j);

        //- Convert a column of a matrix to a Field
        operator Field<cmptType>() const;


    // Member operators

        //- Assignment to a compatible matrix
        template<class Form>
        void operator=(const Matrix<Form, cmptType>&);

        //- Assignment to a compatible const block
        void operator=(const ConstMatrixBlock<MatrixType>&);

        //- Assignment to a compatible block
        void operator=(const MatrixBlock<MatrixType>&);

        //- Assignment to a compatible const block
        template<class MatrixType2>
        void operator=(const ConstMatrixBlock<MatrixType2>&);

        //- Assignment to a compatible block
        template<class MatrixType2>
        void operator=(const MatrixBlock<MatrixType2>&);

        //- Assignment to a compatible MatrixSpace
        template<class MSForm, direction Nrows, direction Ncols>
        void operator=(const MatrixSpace<MSForm, cmptType, Nrows, Ncols>&);

        //- Assignment to a compatible MatrixSpace block
        template
        <
            template<class, direction, direction> class Block,
            class SubTensor,
            direction BRowStart,
            direction BColStart
        >
        void operator=(const Block<SubTensor, BRowStart, BColStart>&);

        //- Assignment to a compatible VectorSpace (column-vector)
        template<class VSForm, direction Ncmpts>
        void operator=(const VectorSpace<VSForm, cmptType, Ncmpts>&);

        //- Assignment to a compatible VectorSpace (column-vector) block
        template
        <
            template<class, direction> class Block,
            class SubVector,
            direction BStart
        >
        void operator=(const Block<SubVector, BStart>&);

        //- Assignment to a Field (column-vector)
        void operator=(const Field<cmptType>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "matrices/MatrixBlock/MatrixBlockI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/MatrixBlock/MatrixBlock.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
