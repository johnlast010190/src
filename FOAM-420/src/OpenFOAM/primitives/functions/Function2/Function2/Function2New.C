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
    (c) 2020 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function2/Constant/Constant2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function2<Type>> Foam::Function2<Type>::New
(
    const word& name,
    const dictionary& dict
)
{
    if (dict.isDict(name))
    {
        const dictionary& coeffsDict(dict.subDict(name));

        const word Function2Type(coeffsDict.lookup("type"));

        const auto ctor = ctorTableLookup("Function2 type", dictionaryConstructorTable_(), Function2Type);
        return ctor(name, coeffsDict);
    }
    else
    {
        Istream& is(dict.lookup(name, false));

        token firstToken(is);
        word Function2Type;

        if (!firstToken.isWord())
        {
            is.putBack(firstToken);
            return autoPtr<Function2<Type>>
            (
                new Function2s::Constant<Type>(name, is)
            );
        }
        else
        {
            Function2Type = firstToken.wordToken();
        }

        const auto ctor = ctorTableLookup("Function2 type", dictionaryConstructorTable_(), Function2Type);
        return ctor(name, dict);
    }
}


// ************************************************************************* //
