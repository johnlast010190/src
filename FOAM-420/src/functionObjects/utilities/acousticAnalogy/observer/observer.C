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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "acousticAnalogy/observer/observer.H"

#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/strings/word/word.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(observer, 0);
    defineRunTimeSelectionTable(observer, word);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::observer>
Foam::observer::New(
                        const word& name,
                        const fvMesh& mesh,
                        const dictionary& dict
                   )
{
    word typeName(dict.lookup("type"));

    const auto ctor = ctorTableLookup("observer type", wordConstructorTable_(), typeName);
    return autoPtr<observer>(ctor(name, mesh, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::observer::observer(
                            const word& name,
                            const fvMesh& mesh,
                            const dictionary& dict
                        )
:
mesh_(mesh),
name_(name),
positions_(),
pPrime_()
{
}


// ************************************************************************* //
