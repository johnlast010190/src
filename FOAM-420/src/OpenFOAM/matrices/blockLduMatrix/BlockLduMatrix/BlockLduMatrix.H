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
    (c) 2019 Esi Ltd.

Class
    BlockLduMatrix

Description
    BlockLduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.  Addressing object must
    be supplied for the upper and lower triangles.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockLduMatrix.C
    BlockLduMatrixOperations.C
    BlockLduMatrixUpdateInterfaces.C
    BlockLduMatrixATmul.C
    BlockLduMatrixHOps.C
    BlockLduMatrixDecouple.C
    BlockLduMatrixDecoupledHOps.C

\*---------------------------------------------------------------------------*/

#ifndef BlockLduMatrix_H
#define BlockLduMatrix_H

#include "fields/CoeffField/coeffFields.H"
#include "fields/FieldFields/FieldField/FieldField.H"
#include "meshes/lduMesh/lduMesh.H"
// #include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduInterface/BlockLduInterfaceField.H"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduInterface/BlockLduInterfaceFieldPtrsList.H"
#include "containers/HashTables/Map/Map.H"
#include "db/IOstreams/IOstreams/InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockLduMatrix;

template<class Type>
Ostream& operator<<(Ostream&, const BlockLduMatrix<Type>&);

template<class Type>
Ostream& operator<<(Ostream&,  const InfoProxy<BlockLduMatrix<Type>>&);

template<class Type>
class BlockConstraint;


/*---------------------------------------------------------------------------*\
                        Class BlockLduMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockLduMatrix
:
    public refCount
{
public:

    // Public data types

        typedef CoeffField<Type> TypeCoeffField;
        typedef Field<Type> TypeField;
        typedef BlockConstraint<Type> ConstraintType;


private:

    // Private data

        // LDU mesh reference
        const lduMesh& lduMesh_;


        // Block matrix elements

            //- Diagonal coefficients
            CoeffField<Type>* diagPtr_;

            //- Upper triangle coefficients.  Also used for symmetric matrix
            CoeffField<Type>* upperPtr_;

            //- Lower triangle coefficients
            CoeffField<Type> *lowerPtr_;


        // Coupling

            //- List of coupled interfaces
            typename BlockLduInterfaceFieldPtrsList<Type>::Type interfaces_;

            //- Coupled interface coefficients, upper
            FieldField<CoeffField, Type> coupleUpper_;

            //- Coupled interface coefficients, lower
            FieldField<CoeffField, Type> coupleLower_;

        // Constraints

            //- Equation triangle map
            mutable Map<ConstraintType> fixedEqns_;


    // Private static data

        //- Matrix constraint fill-in
        //  Equals to the estimated fraction of fixed nodes in the matrix
        static const label fixFillIn;


    // Private member functions

        // Decoupled versions of nmatrix operations

            //- Sum off-diagonal coefficients and add to diagonal,
            //  decoupled version
            void decoupledSumDiag();

            //- Sum negative off-diagonal coefficients and add to diagonal,
            //  decoupled version
            void decoupledNegSumDiag();

            //- Check matrix for diagonal dominance, decoupled version
            void decoupledCheck() const;

            //- Relax matrix, decoupled version
            void decoupledRelax
            (
                const TypeField& x,
                TypeField& b,
                const scalar alpha
            );

            //- Matrix scaling with scalar field, decoupled version
            void decoupledMultEqOp(const scalarField& sf);

            //- Matrix multiplication without coupled interface update,
            //  Decoupled version
            void decoupledAmulCore
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix transpose multiplication without coupled
            //  interface update
            //  Decoupled version
            void decoupledTmulCore
            (
                TypeField& Tx,
                const TypeField& x
            ) const;

            //- Return L-U vector-matrix multiplication in row-form,
            //  Decoupled version
            tmp<TypeField> decoupledH(const TypeField& x) const;

            //- Return L-U  vector-matrix multiplication in off-diagonal form,
            //  Decoupled version
            tmp<TypeField> decoupledFaceH(const TypeField& x) const;


protected:

    // Access to constraints

        //- Return constraint map
        const Map<ConstraintType>& fixedEqns() const
        {
            return fixedEqns_;
        }


        //- Return access constraint map
        Map<ConstraintType>& fixedEqns()
        {
            return fixedEqns_;
        }


public:

    //- Runtime type information
    TypeName("BlockLduMatrix");


    // Constructors

        //- Construct given addressing
        explicit BlockLduMatrix(const lduMesh&);

        //- Construct as copy
        BlockLduMatrix(const BlockLduMatrix<Type>&);

        //- Construct as copy or re-use as specified.
        BlockLduMatrix(BlockLduMatrix<Type>&, bool reUse);

        //- Construct given an LDU addressed mesh and an Istream
        //  from which the coefficients are read
        BlockLduMatrix(const lduMesh&, Istream&);


    // Destructor

        virtual ~BlockLduMatrix();


    // Member functions

        // Access to addressing

            //- Return the LDU mesh from which the addressing is obtained
            virtual const lduMesh& mesh() const
            {
                return lduMesh_;
            }

            //- Return the LDU addressing
            virtual const lduAddressing& lduAddr() const
            {
                return lduMesh_.lduAddr();
            }

            //- Return the patch evaluation schedule
            const lduSchedule& patchSchedule() const
            {
                return lduAddr().patchSchedule();
            }


        // Access to coefficients

            //- Return access to diagonal coefficients
            TypeCoeffField& diag();

            //- Return diagonal coefficients
            const TypeCoeffField& diag() const;

            //- Return access to upper coefficients
            //  Also used for symmetric matrices
            TypeCoeffField& upper();

            //- Return upper coefficients
            //  Also used for symmetric matrices
            const TypeCoeffField& upper() const;

            //- Return access to lower coefficients
            TypeCoeffField& lower();

            //- Return lower coefficients
            const TypeCoeffField& lower() const;

            //- Return access to coupled interface coefficients, upper
            FieldField<CoeffField, Type>& coupleUpper()
            {
                return coupleUpper_;
            }

            //- Return coupled interface coefficients, upper
            const FieldField<CoeffField, Type>& coupleUpper() const
            {
                return coupleUpper_;
            }

            //- Return access to coupled interface coefficients, lower
            FieldField<CoeffField, Type>& coupleLower()
            {
                return coupleLower_;
            }

            //- Return coupled interface coefficients, lower
            const FieldField<CoeffField, Type>& coupleLower() const
            {
                return coupleLower_;
            }

            //- Access to coupled interfaces
            typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces()
            {
                return interfaces_;
            }

            void clearInterfaces();

        // Matrix structure

            //- Return true if there is a diagonal
            bool thereIsDiag() const
            {
                return (diagPtr_);
            }

            //- Return true if upper triangle is allocated
            bool thereIsUpper() const
            {
                return (upperPtr_);
            }

            //- Return true if lower triangle is allocated
            bool thereIsLower() const
            {
                return (lowerPtr_);
            }

            //- Return true if matrix is empty
            bool empty() const;

            //- Return true if matrix is diagonal-only
            bool diagonal() const;

            //- Return true if matrix is symmetric
            bool symmetric() const;

            //- Return true if matrix is asymmetric
            bool asymmetric() const;

            //- Return true if matrix is component-coupled
            bool componentCoupled() const;


        // Operations

            //- Sum off-diagonal coefficients and add to diagonal
            void sumDiag();

            //- Sum negative off-diagonal coefficients and add to diagonal
            void negSumDiag();

            //- Check matrix for diagonal dominance
            void check() const;

            //- Relax and check dominance. Similar to segregated
            void relaxAndCheckDominance
            (
                const TypeField& x,
                TypeField& b,
                const scalar alpha
            );

            //- Relax matrix. Applied only on diagonal elements of diagonal
            void relax
            (
                const TypeField& x,
                TypeField& b,
                const scalar alpha
            );

            //- Relax matrix. Applied only on diagonal elements of diagonal
            void relaxVec
            (
                const TypeField& x,
                TypeField& b,
                const TypeField& alpha
            );

            //- Relax matrix on specific cells
            void relax
            (
                const TypeField& x,
                TypeField& b,
                const labelList& cells,
                const scalar alpha
            );

            //- Matrix multiplication
            void Amul
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix multiplication without coupled interface update
            void AmulCore
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix transpose multiplication
            void Tmul
            (
                TypeField& Ax,
                const TypeField& x
            ) const;

            //- Matrix transpose multiplication without coupled
           //  interface update
            void TmulCore
            (
                TypeField& Ax,
                const TypeField& x
            ) const;


            //- Return decoupled b
            void segregateB
            (
                TypeField& sMul,
                const TypeField& x
            ) const;

            //- Set the residual field using an IOField on the object registry
            //- if it exists
            void setResidualField
            (
                const Field<Type>& residual,
                const word& fieldName
            ) const;


        // Coupled interface functionality

            //- Initialise the update of coupled interfaces
            //  for Amul operations
            void initInterfaces
            (
                const FieldField<CoeffField, Type>& interfaceCoeffs,
                TypeField& result,
                const TypeField& psi
            ) const;

            //- Update coupled interfaces
            void updateInterfaces
            (
                const FieldField<CoeffField, Type>& interfaceCoeffs,
                TypeField& result,
                const TypeField& psi
            ) const;


        // Constraint manipulation

            //- Set constrained value in a prescribed point
            void setValue
            (
                const label eqnIndex,
                const Type& value
            );


        // Residual calculation

            //- Calculate residual
            tmp<TypeField> residual
            (
                const TypeField& x
            ) const;

            tmp<TypeField> residual
            (
                const TypeField& x,
                const TypeField& b
            ) const;

             const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces() const
             {
                   return interfaces_;
             }

        // H-operations

            //- Return L-U vector-matrix multiplication in row-form
            tmp<TypeField> H(const TypeField&) const;

            //- Return L-U  vector-matrix multiplication in off-diagonal form
            tmp<TypeField> faceH(const TypeField&) const;

            //- Partial faceH operation for reducing operations
            //  Only for scalar and vectors currently. It can be easily extended
            tmp<Field<scalar>> faceH
            (
                const direction sI,
                const Field<scalar>&,
                const TypeField&
            ) const;

    // Member operators

        void operator=(const BlockLduMatrix<Type>&);

        void negate();

        void operator+=(const BlockLduMatrix<Type>&);
        void operator-=(const BlockLduMatrix<Type>&);

        void operator*=(const scalarField&);
        void operator*=(const scalar);


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const BlockLduMatrix<Type>&
        );

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const InfoProxy<BlockLduMatrix<Type>>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrix.C"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrixOperations.C"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrixUpdateInterfaces.C"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrixATmul.C"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrixHOps.C"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrixDecouple.C"
#include "matrices/blockLduMatrix/BlockLduMatrix/BlockLduMatrixDecoupledHOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
