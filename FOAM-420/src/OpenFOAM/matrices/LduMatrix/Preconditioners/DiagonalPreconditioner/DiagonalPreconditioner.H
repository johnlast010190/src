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
    Foam::DiagonalPreconditioner

Description
    Diagonal preconditioner for both symmetric and asymmetric matrices.

    The inverse (reciprocal for scalar) of the diagonal is calculated and
    stored.

SourceFiles
    DiagonalPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef DiagonalPreconditioner_H
#define DiagonalPreconditioner_H

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DiagonalPreconditioner Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class DiagonalPreconditioner
:
    public LduMatrix<Type, DType, LUType>::preconditioner
{
    // Private data

        //- The inverse (reciprocal for scalar) diagonal
        Field<DType> rD;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        DiagonalPreconditioner(const DiagonalPreconditioner&);

        //- Disallow default bitwise assignment
        void operator=(const DiagonalPreconditioner&);


public:

    //- Runtime type information
    TypeName("diagonal");


    // Constructors

        //- Construct from matrix components and preconditioner data dictionary
        DiagonalPreconditioner
        (
            const typename LduMatrix<Type, DType, LUType>::solver& sol,
            const dictionary& preconditionerDict
        );


    // Destructor

        virtual ~DiagonalPreconditioner()
        {}


    // Member Functions

        //- Read and reset the preconditioner parameters from the given
        //  dictionary
        virtual void read(const dictionary& preconditionerDict);

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            Field<Type>& wA,
            const Field<Type>& rA
        ) const;

        //- Return wT the transpose-matrix preconditioned form of
        //  residual rT.
        virtual void preconditionT
        (
            Field<Type>& wT,
            const Field<Type>& rT
        ) const
        {
            return(precondition(wT, rT));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matrices/LduMatrix/Preconditioners/DiagonalPreconditioner/DiagonalPreconditioner.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
