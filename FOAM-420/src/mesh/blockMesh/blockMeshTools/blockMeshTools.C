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

#include "blockMeshTools/blockMeshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMeshTools::read
(
    Istream& is,
    label& val,
    const dictionary& dict
)
{
    token t(is);
    if (t.isLabel())
    {
        val = t.labelToken();
    }
    else if (t.isWord())
    {
        const word& varName = t.wordToken();
        const entry* ePtr = dict.lookupScopedEntryPtr
        (
            varName,
            true,
            true
        );
        if (ePtr)
        {
            // Read as label
            val = Foam::readLabel(ePtr->stream());
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Undefined variable "
                << varName << ". Valid variables are " << dict
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Illegal token " << t.info()
            << " when trying to read label"
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);
}


Foam::label Foam::blockMeshTools::read
(
    Istream& is,
    const dictionary& dict
)
{
    label val;
    read(is, val, dict);
    return val;
}


void Foam::blockMeshTools::write
(
    Ostream& os,
    const label val,
    const dictionary& dict
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isStream())
        {
            label keyVal(Foam::readLabel(iter().stream()));
            if (keyVal == val)
            {
                os << iter().keyword();
                return;
            }
        }
    }
    os << val;
}


const Foam::keyType& Foam::blockMeshTools::findEntry
(
    const dictionary& dict,
    const label val
)
{
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isStream())
        {
            label keyVal(Foam::readLabel(iter().stream()));
            if (keyVal == val)
            {
                return iter().keyword();
            }
        }
    }

    return keyType::null;
}


// ************************************************************************* //
