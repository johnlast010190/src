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

InNamespace
    Foam

Description
    Optimized hashing functions.

See also
    Hasher.H

\*---------------------------------------------------------------------------*/

#ifndef HasherInt_H
#define HasherInt_H

#include "primitives/hashes/Hasher/Hasher.H"
#include "primitives/ints/int/int.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- An optimized version of Hasher
    //  \param[in] data - an array of uint32_t values
    //  \param[in] length - the number of values (not bytes)
    //  \param[in] seed - the previous hash, or an arbitrary value
    unsigned HasherInt(const uint32_t* data, size_t length, unsigned seed = 0);

    //- An optimized version of Hasher, returning dual hash values
    //  \param[in] data - an array of uint32_t values
    //  \param[in] length - the number of values (not bytes)
    //  \param[in] hash1 - the previous hash, or an arbitrary value
    //    on output, the primary hash value
    //  \param[in] hash1 - the previous hash, or an arbitrary value
    //    on output, the secondary hash value
    unsigned HasherDual
    (
        const uint32_t* data,
        size_t length,
        unsigned& hash1,
        unsigned& hash2
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
