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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/Tensor/floatTensor/floatTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if !defined(FOAM_SP)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char* const Foam::floatTensor::vsType::typeName = "floatTensor";

template<>
const char* const Foam::floatTensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
    "yx", "yy", "yz",
    "zx", "zy", "zz"
};

template<>
const Foam::floatTensor Foam::floatTensor::vsType::zero
(
    floatTensor::uniform(0)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::one
(
    floatTensor::uniform(1)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::max
(
    floatTensor::uniform(floatScalarVGREAT)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::min
(
    floatTensor::uniform(-floatScalarVGREAT)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::rootMax
(
    floatTensor::uniform(floatScalarROOTVGREAT)
);

template<>
const Foam::floatTensor Foam::floatTensor::vsType::rootMin
(
    floatTensor::uniform(-floatScalarROOTVGREAT)
);

template<>
const Foam::floatTensor Foam::floatTensor::I
(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif //!defined(FOAM_SP)

// ************************************************************************* //
