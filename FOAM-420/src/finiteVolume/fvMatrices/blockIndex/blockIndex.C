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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMatrices/blockIndex/blockIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label blockIndex::index
(
    const word& nameField
) const
{
    label size = names_.size();

    int i = 0;
    bool found = false;
    label index = -1;

    do
    {
        if (names_[i] == nameField)
        {
            index = i;
            found = true;
        }
        i += 1;

    } while (!found && i < size);

    return index;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

blockIndex::blockIndex
(
    const wordList& names,
    const labelList& sizes
)
:
    names_(names),
    sizes_(sizes),
    startIndices_(names.size(), 0),
    gNames_(0)
{
    if (names.size() > 0)
    {
        startIndices_[0] = 0;
        for (label i = 1; i < names.size(); i++)
        {
            startIndices_[i] = startIndices_[i-1] + sizes_[i-1];
        }
    }
    forAll(names, nI)
    {
        for (direction sI = 0; sI < sizes_[nI]; sI++)
        {
            gNames_.append(names[nI]);
        }
    }
}


blockIndex::blockIndex
(
    const blockIndex& bI
)
:
    names_(bI.names_),
    sizes_(bI.sizes_),
    startIndices_(bI.startIndices_),
    gNames_(bI.gNames_)
{}


blockIndex::blockIndex()
:
    names_(0),
    sizes_(0),
    startIndices_(0),
    gNames_(0)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void blockIndex::check()
{
    if
    (
        names_.size() == 0
     || sizes_.size() == 0
     || startIndices_.size() == 0
    )
    {
        FatalErrorInFunction
            << "blockIndex is not setup for the fvBlockMatrix"
            << abort(FatalError);
    }
}


label blockIndex::findIndex
(
    const word& nameField
) const
{
    return index(nameField);
}


label blockIndex::findStartIndex
(
    const word& nameField
) const
{
    return startIndices_[index(nameField)];
}


bool blockIndex::increment
(
    const label i1,
    const label i2
) const
{
    return startIndices_[i1]<startIndices_[i2];
}


word blockIndex::variableName
(
    const label index
) const
{
    return gNames_[index];
}


void blockIndex::operator=
(
    const blockIndex& bI
)
{
    names_ = bI.names_;
    sizes_ = bI.sizes_;
    startIndices_ = bI.startIndices_;
    gNames_ = bI.gNames_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
