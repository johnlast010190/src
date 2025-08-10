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

\*---------------------------------------------------------------------------*/

#include "matrices/blockLduMatrix/BlockLduPrecons/BlockLduPrecon/BlockLduPrecon.H"
#include "matrices/blockLduMatrix/BlockLduPrecons/BlockNoPrecon/blockNoPrecons.H"

template<class Type>
class BlockNoPrecon;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduPrecon<Type>> Foam::BlockLduPrecon<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const word fieldName
)
{
    word preconName = getName(dict);

    const dictionary& controls = dict.subOrEmptyDict(preconName+"Coeffs");

    if (matrix.diagonal())
    {
        // No preconditioning for the diagonal matrix
        return autoPtr<BlockLduPrecon<Type>>
        (
            new BlockNoPrecon<Type>
            (
                matrix,
                controls,
                fieldName
            )
        );
    }
    else
    {
        typename dictionaryConstructorTable::iterator constructorIter =
            dictionaryConstructorTable_().find(preconName);

        if (constructorIter == dictionaryConstructorTable_().end())
        {
            FatalIOErrorIn
            (
                "autoPtr<BlockLduPrecon> BlockLduPrecon::New\n"
                "(\n"
                "    const BlockLduMatrix<Type>& matrix,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown matrix preconditioner " << preconName
                << endl << endl
                << exit(FatalIOError);
        }

        return autoPtr<BlockLduPrecon<Type>>
        (
            constructorIter->second
            (
                matrix,
                controls,
                fieldName
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::BlockLduPrecon<Type>::getName(const dictionary& dict)
{
    word name;

    dict.lookup("preconditioner") >> name;

    return name;
}


// ************************************************************************* //
