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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2014 blueCAPE Lda

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.

    Modifications made:
      - Derived from the patches for blueCFD 2.1 and 2.2.
      - The method 'Random::bit()' now relies on the system dependent method
        'osRandomBit()'.

\*---------------------------------------------------------------------------*/

#include "primitives/random/Random/Random.H"
#include "include/OSspecific.H"
#include "db/IOstreams/Pstreams/PstreamReduceOps.H"

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * //

#if !defined(WIN64) && !defined(WIN32)
Foam::scalar Foam::Random::scalar01()
{
    return osRandomDouble(buffer_);
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Random::Random(const label seed)
:
    #if !defined(WIN64) && !defined(WIN32)
    buffer_(osRandomBufferSize()),
    #endif
    seed_(seed),
    hasGaussSample_(false),
    gaussSample_(0)
{
    // Initialise the random number generator
#if !defined(WIN64) && !defined(WIN32)
    osRandomSeed(seed_, buffer_);
    #endif
}


Foam::Random::Random(const Random& r, const bool reset)
:
#if !defined(WIN64) && !defined(WIN32)
    buffer_(r.buffer_),
    #endif
    seed_(r.seed_),
    hasGaussSample_(r.hasGaussSample_),
    gaussSample_(r.gaussSample_)
{
    if (reset)
    {
        hasGaussSample_ = false;
        gaussSample_ = 0;

        // Re-initialise the samples
 #if !defined(WIN64) && !defined(WIN32)
        osRandomSeed(seed_, buffer_);
        #else
        osRandomSeed(seed_);
        #endif
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Random::~Random()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Random::reset(const label seed)
{
    seed_ = seed;
#if !defined(WIN64) && !defined(WIN32)
    osRandomSeed(seed_, buffer_);
    #else
    osRandomSeed(seed_);
    #endif
}


int Foam::Random::bit()
{
#if !defined(WIN64) && !defined(WIN32)
    return osRandomInteger(buffer_) & 1;
    #else
    return osRandomInteger() & 1;
}

Foam::scalar Foam::Random::scalar01()
{
    return osRandomDouble();
}


Foam::vector Foam::Random::vector01()
{
    Foam::vector rndVec;
    for (direction cmpt=0; cmpt<Foam::vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) = scalar01();
    }

    return rndVec;
    #endif
}


template<>
Foam::scalar Foam::Random::sample01()
{
    return scalar01();
}


template<>
Foam::label Foam::Random::sample01()
{
    return round(scalar01());
}


template<>
Foam::scalar Foam::Random::GaussNormal()
{
    if (hasGaussSample_)
    {
        hasGaussSample_ = false;
        return gaussSample_;
    }
    else
    {
        scalar rsq, v1, v2;
        do
        {
            v1 = 2*scalar01() - 1;
            v2 = 2*scalar01() - 1;
            rsq = sqr(v1) + sqr(v2);
        } while (rsq >= 1 || rsq == 0);

        scalar fac = sqrt(-2*log(rsq)/rsq);

        gaussSample_ = v1*fac;
        hasGaussSample_ = true;

        return v2*fac;
    }
}


template<>
Foam::label Foam::Random::GaussNormal()
{
    return round(GaussNormal<scalar>());
}


template<>
Foam::scalar Foam::Random::position
(
    const scalar& start,
    const scalar& end
)
{
    return start + scalar01()*(end - start);
}


template<>
Foam::label Foam::Random::position(const label& start, const label& end)
{
    return start + round(scalar01()*(end - start));
}


template<>
Foam::scalar Foam::Random::globalSample01()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01();
    }

    Pstream::scatter(value);

    return value;
}


template<>
Foam::label Foam::Random::globalSample01()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01();
    }

    Pstream::scatter(value);

    return round(value);
}


template<>
Foam::scalar Foam::Random::globalGaussNormal()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = GaussNormal<scalar>();
    }

    Pstream::scatter(value);

    return value;
}


template<>
Foam::label Foam::Random::globalGaussNormal()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = GaussNormal<scalar>();
    }

    Pstream::scatter(value);

    return round(value);
}


template<>
Foam::scalar Foam::Random::globalPosition
(
    const scalar& start,
    const scalar& end
)
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01()*(end - start);
    }

    Pstream::scatter(value);

    return start + value;
}


template<>
Foam::label Foam::Random::globalPosition
(
    const label& start,
    const label& end
)
{
    label value = labelMin;

    if (Pstream::master())
    {
        value = round(scalar01()*(end - start));
    }

    Pstream::scatter(value);

    return start + value;
}


// ************************************************************************* //
