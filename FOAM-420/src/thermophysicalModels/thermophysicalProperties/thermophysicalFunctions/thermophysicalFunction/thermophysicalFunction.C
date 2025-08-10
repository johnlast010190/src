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

\*---------------------------------------------------------------------------*/

#include "thermophysicalFunctions/thermophysicalFunction/thermophysicalFunction.H"
#include "containers/HashTables/HashTable/HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalFunction, 0);
    defineRunTimeSelectionTable(thermophysicalFunction, Istream);
    defineRunTimeSelectionTable(thermophysicalFunction, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermophysicalFunction> Foam::thermophysicalFunction::New
(
    Istream& is
)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing thermophysicalFunction"
            << endl;
    }

    const word thermophysicalFunctionType(is);

    const auto ctor = ctorTableLookup("thermophysicalFunction type", IstreamConstructorTable_(), thermophysicalFunctionType);
    return autoPtr<thermophysicalFunction>(ctor(is));
}


Foam::autoPtr<Foam::thermophysicalFunction> Foam::thermophysicalFunction::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction
            << "Constructing thermophysicalFunction"
            << endl;
    }

    const word thermophysicalFunctionType(dict.lookup("functionType"));

    const auto ctor = ctorTableLookup("thermophysicalFunction type", dictionaryConstructorTable_(), thermophysicalFunctionType);
    return autoPtr<thermophysicalFunction>(ctor(dict));
}


// ************************************************************************* //
