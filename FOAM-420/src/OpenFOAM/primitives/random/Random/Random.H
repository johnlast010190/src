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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2020 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::Random

Description
    Random number generator.


SourceFiles
    RandomI.H
    Random.C
    RandomTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Random_H
#define Random_H

#include "containers/Lists/List/List.H"
#if defined(WIN64) || defined(WIN32)
#include "include/OSspecific.H"
#include "primitives/Vector/vector/vector.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Random Declaration
\*---------------------------------------------------------------------------*/

class Random
{
    // Private data

        //- Buffer used by re-entrant random number generator
#if !defined(WIN64) && !defined(WIN32)
        List<char> buffer_;
#endif

        //- Initial random number seed
        label seed_;

        //- Indicator, which tells if there is a stored gaussian sample
        bool hasGaussSample_;

        //- Stored sample value
        scalar gaussSample_;


    // Private Member Functions

#if !defined(WIN64) && !defined(WIN32)
        //- Returns the current sample
        scalar scalar01();
#endif


public:

    // Constructors

        //- Construct given seed and sample count
        explicit Random(const label seed = 123456);

        //- Copy constructor with optional reset of sampleI
        Random(const Random& r, const bool reset);


    // Destructor
    ~Random();


    // Member functions

        // Access

            //- Return const access to the initial random number seed
            inline label seed() const;

            //- Reset the random number generator seed
            void reset(const label seed);


        // Evaluation

            // Random numbers

                //- Return a random bit
                int bit();
#if defined(WIN64) || defined(WIN32)
                //- Scalar [0..1] (so including 0,1)
                Foam::scalar scalar01();

                //- Vector with every component scalar01
                Foam::vector vector01();
#endif

                //- Return a sample whose components lie in the range 0-1
                template<class Type>
                Type sample01();

                //- Return a sample whose components are normally distributed
                //  with zero mean and unity variance N(0, 1)
                template<class Type>
                Type GaussNormal();

                //- Return a sample between start and end
                template<class Type>
                Type position(const Type& start, const Type& end);

                //- Randomise value in the range 0-1
                template<class Type>
                void randomise01(Type& value);

                //- Shuffle the values in the list
                template<class Type>
                void shuffle(UList<Type>& values);


            // Global random numbers - consistent across all processors

                //- Return a sample whose components lie in the range 0-1
                template<class Type>
                Type globalSample01();

                //- Return a sample whose components are normally distributed
                //  with zero mean and unity variance N(0, 1)
                template<class Type>
                Type globalGaussNormal();

                //- Return a sample between start and end
                template<class Type>
                Type globalPosition(const Type& start, const Type& end);

                //- Randomise value in the range 0-1
                template<class Type>
                void globalRandomise01(Type& value);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Template specialisations

template<>
scalar Random::sample01<scalar>();

template<>
label Random::sample01<label>();

template<>
scalar Random::GaussNormal<scalar>();

template<>
label Random::GaussNormal<label>();

template<>
scalar Random::position<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
label Random::position<label>(const label& start, const label& end);

template<>
scalar Random::globalSample01<scalar>();

template<>
label Random::globalSample01<label>();

template<>
scalar Random::globalGaussNormal<scalar>();

template<>
label Random::globalGaussNormal<label>();

template<>
scalar Random::globalPosition<scalar>
(
    const scalar& start,
    const scalar& end
);

template<>
label Random::globalPosition<label>(const label& start, const label& end);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "primitives/random/Random/RandomI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "primitives/random/Random/RandomTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
