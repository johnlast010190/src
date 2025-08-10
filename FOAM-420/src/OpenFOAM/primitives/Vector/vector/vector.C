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

Description
    Vector of scalars.

\*---------------------------------------------------------------------------*/

#include "primitives/Vector/vector/vector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::vector::vsType::typeName = "vector";

template<>
const char* const Foam::vector::vsType::componentNames[] = {"x", "y", "z"};

template<>
const Foam::vector Foam::vector::vsType::zero(vector::uniform(0));

template<>
const Foam::vector Foam::vector::vsType::one(vector::uniform(1));

template<>
const Foam::vector Foam::vector::vsType::max(vector::uniform(VGREAT));

template<>
const Foam::vector Foam::vector::vsType::min(vector::uniform(-VGREAT));

template<>
const Foam::vector Foam::vector::vsType::rootMax(vector::uniform(ROOTVGREAT));

template<>
const Foam::vector Foam::vector::vsType::rootMin(vector::uniform(-ROOTVGREAT));


// ************************************************************************* //
