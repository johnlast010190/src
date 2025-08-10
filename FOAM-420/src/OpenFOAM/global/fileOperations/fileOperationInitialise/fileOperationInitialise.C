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
    (c) 2017-2018 OpenFOAM Foundation
    (c) 2019 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "global/fileOperations/fileOperationInitialise/fileOperationInitialise.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(fileOperationInitialise, 0);
    defineRunTimeSelectionTable(fileOperationInitialise, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::fileOperationInitialise::fileOperationInitialise
(
    int& argc,
    char**& argv
)
{}


Foam::autoPtr<Foam::fileOperations::fileOperationInitialise>
Foam::fileOperations::fileOperationInitialise::New
(
    const word& type,
    int& argc,
    char**& argv
)
{
    DebugInFunction << "Constructing fileOperationInitialise" << endl;

    const auto ctor = ctorTableLookup("fileOperationInitialise type", wordConstructorTable_(), type);
    return autoPtr<fileOperationInitialise>(ctor(argc, argv));
}


// ************************************************************************* //
