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

\*---------------------------------------------------------------------------*/

#include "primitives/random/random11/random11.H"
#include <initializer_list>
#if defined(WIN64) || defined(WIN32)
#include "include/OSspecific.H"
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::random11::random11(int seed)
:
    randomEngine_(seed)
{}

// * * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

int Foam::random11::bit()
{
    std::uniform_int_distribution<int> dist(0, 1);
    return dist(randomEngine_);
}


Foam::scalar Foam::random11::scalar01()
{
    std::uniform_real_distribution<scalar> dist(0, 1);
    return dist(randomEngine_);
}


Foam::vector Foam::random11::vector01()
{
    vector rndVec;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) = scalar01();
    }

    return rndVec;
}


Foam::sphericalTensor Foam::random11::sphericalTensor01()
{
    sphericalTensor rndTen;
    rndTen.ii() = scalar01();

    return rndTen;
}


Foam::symmTensor Foam::random11::symmTensor01()
{
    symmTensor rndTen;
    for (direction cmpt=0; cmpt<symmTensor::nComponents; cmpt++)
    {
        rndTen.component(cmpt) = scalar01();
    }

    return rndTen;
}


Foam::tensor Foam::random11::tensor01()
{
    tensor rndTen;
    for (direction cmpt=0; cmpt<tensor::nComponents; cmpt++)
    {
        rndTen.component(cmpt) = scalar01();
    }

    return rndTen;
}


Foam::label Foam::random11::integer(const label lower, const label upper)
{
    std::uniform_int_distribution<label> dist(lower, upper);
    return dist(randomEngine_);
}


Foam::vector Foam::random11::position(const vector& start, const vector& end)
{
    vector rndVec(start);

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) +=
            scalar01()*(end.component(cmpt) - start.component(cmpt));
    }

    return rndVec;
}


void Foam::random11::randomise(scalar& s)
{
     s = scalar01();
}


void Foam::random11::randomise(vector& v)
{
    v = vector01();
}


void Foam::random11::randomise(sphericalTensor& st)
{
    st = sphericalTensor01();
}


void Foam::random11::randomise(symmTensor& st)
{
    st = symmTensor01();
}


void Foam::random11::randomise(tensor& t)
{
    t = tensor01();
}


Foam::scalar Foam::random11::GaussNormal()
{
    std::normal_distribution<scalar> dist;
    return dist(randomEngine_);
}

Foam::scalar Foam::random11::randomise(scalar a, scalar b)
{
    std::uniform_real_distribution<scalar> dist(a, b);
    return dist(randomEngine_);
}

Foam::label Foam::random11::randomise(label a, label b)
{
    std::uniform_int_distribution<label> dist(a, b);
    return dist(randomEngine_);
}

Foam::label Foam::random11::randomise(const List<scalar>& weights)
{
    std::discrete_distribution<label> dist(weights.begin(), weights.end());
    return dist(randomEngine_);
}

// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
