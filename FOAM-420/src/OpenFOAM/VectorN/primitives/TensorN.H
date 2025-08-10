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
    (c) 2016 Esi Ltd.
    (c) 2010 Ivor Clifford

Class
    TensorN

Description

SourceFiles
    TensorNI.H

\*---------------------------------------------------------------------------*/

#ifndef TensorN_H
#define TensorN_H

#include "VectorN/primitives/VectorN.H"
#include "VectorN/primitives/SphericalTensorN.H"
#include "VectorN/primitives/DiagTensorN.H"
#include "containers/Lists/List/ListDecl.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

template <class Cmpt, int length>
class DiagTensorN;

template <class Cmpt, int length>
class SphericalTensorN;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class TensorN Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt, int length>
class TensorN
:
    public VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>
{

public:

    // Member constants

        constexpr static int rank = 2;
        constexpr static int rowLength = length;

    // Static data members

        static const char* const typeName;

        static const TensorN zero;
        static const TensorN one;


    // Constructors

        //- Construct null
        inline TensorN();

        //- Construct initialized to zero
        inline TensorN(const Foam::zero);

        //- Construct given VectorSpace
        inline TensorN
        (
            const VectorSpace<TensorN<Cmpt, length>, Cmpt, length*length>&
        );

        //- Construct given component value.  Special use only!
        explicit inline TensorN(const Cmpt& tx);

        //- Construct from Istream
        TensorN(Istream&);


    // Member Functions

        //- Return (i, j) component
        inline const Cmpt& operator()
        (
            const direction i,
            const direction j
        ) const;

        //- Return access to (i, j) component
        inline Cmpt& operator()
        (
            const direction i,
            const direction j
        );

        //- Diagonal
        inline DiagTensorN<Cmpt, length> diag() const;

        //- Transpose
        inline TensorN<Cmpt, length> T() const;

        //- Negative sum the vertical off-diagonal components
        inline TensorN<Cmpt, length> negSumDiag() const;

    // Member Operators

        //- Assign to a zero
        inline void operator=(const Foam::zero&);

        //- Assign to a SphericalTensorN
        inline void operator=(const SphericalTensorN<Cmpt, length>&);

        //- Assign to a DiagTensorN
        inline void operator=(const DiagTensorN<Cmpt, length>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "VectorN/primitives/TensorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
