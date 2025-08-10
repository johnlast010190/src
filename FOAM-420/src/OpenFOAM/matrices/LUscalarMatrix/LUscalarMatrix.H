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

Class
    Foam::LUscalarMatrix

Description
    Class to perform the LU decomposition on a symmetric matrix.

SourceFiles
    LUscalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef LUscalarMatrix_H
#define LUscalarMatrix_H

#include "matrices/scalarMatrices/scalarMatrices.H"
#include "primitives/ints/lists/labelList.H"
#include "fields/FieldFields/FieldField/FieldField.H"
#include "matrices/lduMatrix/lduAddressing/lduInterfaceFields/lduInterfaceField/lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMatrix;
class procLduMatrix;

/*---------------------------------------------------------------------------*\
                           Class LUscalarMatrix Declaration
\*---------------------------------------------------------------------------*/

class LUscalarMatrix
:
    public scalarSquareMatrix
{
    // Private data

        //- Communicator to use
        const label comm_;

        //- Processor matrix offsets
        labelList procOffsets_;

        //- The pivot indices used in the LU decomposition
        labelList pivotIndices_;


    // Private member functions

        //- Convert the given lduMatrix into this LUscalarMatrix
        void convert
        (
            const lduMatrix& ldum,
            const FieldField<Field, scalar>& interfaceCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );

        //- Convert the given list of procLduMatrix into this LUscalarMatrix
        //  on the master processor
        void convert(const PtrList<procLduMatrix>& lduMatrices);


        //- Print the ratio of the mag-sum of the off-diagonal coefficients
        //  to the mag-diagonal
        void printDiagonalDominance() const;


public:

    // Declare name of the class and its debug switch
    ClassName("LUscalarMatrix");


    // Constructors

        //- Construct null
        LUscalarMatrix();

        //- Construct from and perform LU decomposition of the matrix M
        LUscalarMatrix(const scalarSquareMatrix& M);

        //- Construct from lduMatrix and perform LU decomposition
        LUscalarMatrix
        (
            const lduMatrix&,
            const FieldField<Field, scalar>& interfaceCoeffs,
            const lduInterfaceFieldPtrsList& interfaces
        );


    // Member Functions

        //- Perform the LU decomposition of the matrix M
        void decompose(const scalarSquareMatrix& M);

        //- Solve the linear system with the given source
        //  and returning the solution in the Field argument x.
        //  This function may be called with the same field for x and source.
        template<class Type>
        void solve(Field<Type>& x, const Field<Type>& source) const;

        //- Solve the linear system with the given source
        //  returning the solution
        template<class Type>
        tmp<Field<Type>> solve(const Field<Type>& source) const;

        //- Set M to the inverse of this square matrix
        void inv(scalarSquareMatrix& M) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/LUscalarMatrix/LUscalarMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
