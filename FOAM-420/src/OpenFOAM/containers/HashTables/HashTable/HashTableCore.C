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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/ints/uLabel/uLabel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(HashTableCore, 0);
}

// Approximately labelMax/4
const Foam::label Foam::HashTableCore::maxTableSize(1L << (sizeof(label)*8-3));

const Foam::HashTableCore::NullIter Foam::HashTableCore::NullIter::nullIter;
const Foam::HashTableCore::NullIter*
Foam::HashTableCore::NullIter::nullIterPtr
(
    &Foam::HashTableCore::NullIter::nullIter
);


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::HashTableCore::canonicalSize(const label requested_size)
{
    if (requested_size < 1)
    {
        return 0;
    }
    else if (requested_size >= maxTableSize)
    {
        return maxTableSize;
    }

    // Enforce power of two - makes for a very fast modulus.
    // Use unsigned for these calculations.
    //
    // - The lower limit (8) is somewhat arbitrary, but if the hash table
    //   is too small, there will be many direct table collisions.
    // - The upper limit (approx. labelMax/4) must be a power of two,
    //   need not be extremely large for hashing.

    uLabel powerOfTwo = 8u; // lower-limit

    const uLabel size = requested_size;
    if (size <= powerOfTwo)
    {
        return powerOfTwo;
    }
    else if (size & (size-1))  // <- Modulus of i^2
    {
        // Determine power-of-two. Brute-force is fast enough.
        while (powerOfTwo < size)
        {
            powerOfTwo <<= 1;
        }

        return powerOfTwo;
    }
    else
    {
        return size;
    }
}


// ************************************************************************* //
