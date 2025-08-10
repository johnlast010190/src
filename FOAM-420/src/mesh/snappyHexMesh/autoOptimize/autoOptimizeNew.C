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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"
#include "autoOptimize/autoOptimize.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::autoOptimize> Foam::autoOptimize::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    word optimMethod = dict.lookupOrDefault<word>("type", "cfMeshOptimize");

    Info<< "Using optimization method " << optimMethod << endl;

    const auto ctor = ctorTableLookup("autoOptimize type", dictionaryConstructorTable_(), optimMethod);
    return autoPtr<autoOptimize>(ctor(mesh, dict));
}

// ************************************************************************* //
