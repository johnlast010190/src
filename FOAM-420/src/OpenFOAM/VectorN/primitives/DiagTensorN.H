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
    DiagTensorN

Description
    Templated NXN DiagTensor derived from VectorSpace adding construction from
    N components, and the inner-product (dot-product) and outer-product operators.

SourceFiles
    DiagTensorNI.H

\*---------------------------------------------------------------------------*/

#ifndef DiagTensorN_H
#define DiagTensorN_H

#include "VectorN/primitives/VectorN.H"
#include "VectorN/primitives/SphericalTensorN.H"
#include "containers/Lists/List/ListDecl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DiagTensor Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt, int length>
class DiagTensorN
:
    public VectorSpace<DiagTensorN<Cmpt, length>, Cmpt, length>
{

public:

    // Member constants

        constexpr static int rank = 2;
        constexpr static int rowLength = length;


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const DiagTensorN zero;
        static const DiagTensorN one;
        static const DiagTensorN I;


    // Constructors

        //- Construct null
        inline DiagTensorN();

        //- Construct initialized to zero
        inline DiagTensorN(const Foam::zero);

        //- Construct given VectorSpace
        inline DiagTensorN(const VectorSpace<DiagTensorN<Cmpt, length>, Cmpt, length>&);

        //- Construct from Istream
        inline DiagTensorN(Istream&);

        //- Construct given component value.  Special use only!
        explicit inline DiagTensorN(const Cmpt& tx);


    // Member Functions

        //- Diagonal
        inline DiagTensorN<Cmpt, length> diag() const;

        //- Transpose
        inline DiagTensorN<Cmpt, length> T() const;

    // Member Operators

        //- Assign to a zero
        inline void operator=(const Foam::zero&);

        //- Assign to a SphericalTensorN
        inline void operator=(const SphericalTensorN<Cmpt, length>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "VectorN/primitives/DiagTensorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
