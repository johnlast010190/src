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
    (c) 2004-6 H. Jasak All rights reserved

Class
    BlockLduSmoother

Description
    Block LDU matrix smoother virtual base class

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduSmoothers/BlockLduSmoother/BlockLduSmoother.H"

template<class Type>
class BlockNoSmoother;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduSmoother<Type>> Foam::BlockLduSmoother<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
{
    word smootherName;

    // Handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("smoother", false, false);
    if (e.isDict())
    {
        e.dict().lookup("smoother") >> smootherName;
    }
    else
    {
        e.stream() >> smootherName;
    }

    // Not (yet?) needed:
    // const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    typename dictionaryConstructorTable::iterator constructorIter =
        dictionaryConstructorTable_().find(smootherName);

    if (constructorIter == dictionaryConstructorTable_().end())
    {
        FatalIOErrorIn
        (
            "autoPtr<BlockLduSmoother> BlockLduSmoother::New\n"
            "(\n"
            "    const BlockLduMatrix<Type>& matrix,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown matrix smoother " << smootherName
            << endl << endl
            << exit(FatalIOError);
    }

    return autoPtr<BlockLduSmoother<Type>>
    (
        constructorIter->second
        (
            matrix,
            dict
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::BlockLduSmoother<Type>::getName(const dictionary& dict)
{
    word name;

    // handle primitive or dictionary entry
    const entry& e = dict.lookupEntry("preconditioner", false, false);
    if (e.isDict())
    {
        e.dict().lookup("preconditioner") >> name;
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


// ************************************************************************* //
