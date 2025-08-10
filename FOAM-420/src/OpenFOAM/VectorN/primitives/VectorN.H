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
    VectorN

Description
    Templated 2D Vector derived from VectorSpace adding construction from
    2 components, element access using x() and y() member functions and
    the inner-product (dot-product).

SourceFiles
    VectorNI.H

\*---------------------------------------------------------------------------*/

#ifndef VectorN_H
#define VectorN_H

#include "primitives/VectorSpace/VectorSpace.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/zero/zero.H"
#include "containers/Lists/List/ListDecl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class VectorN Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt, int length>
class VectorN
:
    public VectorSpace<VectorN<Cmpt, length>, Cmpt, length>
{

public:

    // Member constants

        enum
        {
            rank = 1 // Rank of VectorN is 1
        };


    // Static data members

        static const char* const typeName;
        static const VectorN zero;
        static const VectorN one;
        static const VectorN max;
        static const VectorN min;


    // Constructors

        //- Construct null
        inline VectorN();

        //- Construct initialized to zero
        inline VectorN(const Foam::zero);

        //- Construct given VectorSpace
        inline VectorN
        (
            const VectorSpace<VectorN<Cmpt, length>, Cmpt, length>&
        );

        //- Construct given component value.  Special use only!
        explicit inline VectorN(const Cmpt& vx);

        //- Construct from Istream
        inline VectorN(Istream&);


    // Member Functions

        //- Return ith component
        inline const Cmpt& operator()(const direction i) const;

        //- Return access to ith component
        inline Cmpt& operator()(const direction i);

        //- Componentwise multiply
        inline VectorN<Cmpt, length> cmptMultiply
        (
            const VectorN<Cmpt, length>&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "VectorN/primitives/VectorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
