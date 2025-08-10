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
    (c) 2017 Esi Ltd.

Class
    Foam::random11

Description
    Random number generator using c++11 with Mersenne Twister engine.

SourceFiles
    random11.C

\*---------------------------------------------------------------------------*/

#ifndef FOAM_RANDOM11_H
#define FOAM_RANDOM11_H

#include <random>
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/Vector/vector/vector.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/SphericalTensor/sphericalTensor/sphericalTensor.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "containers/Lists/List/List.H"

namespace Foam
{

class random11
{
    //- random number generator
    std::mt19937 randomEngine_;

public:

    // Constructors

        //- construct from seed
        random11(int seed=1);

    // Member functions

        int bit();

        //- Scalar [0..1] (so including 0,1)
        scalar scalar01();

        //- Vector with every component scalar01
        vector vector01();

        //- sphericalTensor with every component scalar01
        sphericalTensor sphericalTensor01();

        //- symmTensor with every component scalar01
        symmTensor symmTensor01();

        //- Tensor with every component scalar01
        tensor tensor01();

        //- Label [lower..upper]
        label integer(const label lower, const label upper);

        vector position(const vector&, const vector&);

        void randomise(scalar&);
        void randomise(vector&);
        void randomise(sphericalTensor&);
        void randomise(symmTensor&);
        void randomise(tensor&);

        // range randomisers

            //-random number between specified scalars
            scalar randomise(scalar, scalar);

            //-random number between specified label
            label randomise(label, label);

            // weighted random integer between 0 - n (list size)
            label randomise(const List<scalar>& weights);

        //- Return a normal Gaussian randon number
        //  with zero mean and unity variance N(0, 1)
        scalar GaussNormal();
};

} //end namespace foam

#endif // FOAM_RANDOM11_H
