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
    (c) 2019 Esi Ltd.

Class
    Foam::SphericalTensor

Description
    Templated 3D SphericalTensor derived from VectorSpace adding construction
    from 1 component, element access using th ii() member function and the
    inner-product (dot-product) and outer-product operators.

SourceFiles
    SphericalTensorI.H

\*---------------------------------------------------------------------------*/

#ifndef SphericalTensor_H
#define SphericalTensor_H

#include "primitives/VectorSpace/VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SphericalTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SphericalTensor
:
    public VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef SphericalTensor<label> labelType;


    // Member constants

        //- Rank of SphericalTensor is 2
        static const direction rank = 2;


    // Static data members

        static const SphericalTensor I;
        static const SphericalTensor oneThirdI;
        static const SphericalTensor twoThirdsI;


    //- Component labeling enumeration
    enum components { II };


    // Constructors

        //- Construct null
        inline SphericalTensor();

        //- Construct initialized to zero
        inline SphericalTensor(const Foam::zero);

        //- Construct given VectorSpace
        template<class Cmpt2>
        inline SphericalTensor
        (
            const VectorSpace<SphericalTensor<Cmpt2>, Cmpt2, 1>&
        );

        //- Construct given the component
        inline SphericalTensor(const Cmpt& tii);

        //- Construct from Istream
        inline SphericalTensor(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& ii() const;
            inline Cmpt& ii();


        //- Transpose
        inline const SphericalTensor<Cmpt>& T() const;

        //- Operators
        inline const Cmpt& operator()(const direction&, const direction&) const;

        //-
        inline Cmpt& operator()(const direction&, const direction&);

        //-
        inline const Cmpt& operator()(const direction) const;

        //-
        inline Cmpt& operator()(const direction);

};

template<typename Cmpt>
constexpr bool is_contiguous<SphericalTensor<Cmpt>> = contiguous<VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>>();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "primitives/SphericalTensor/SphericalTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
