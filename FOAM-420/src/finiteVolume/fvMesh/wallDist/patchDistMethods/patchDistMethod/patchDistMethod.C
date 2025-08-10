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
    (c) 2017 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fvMesh/wallDist/patchDistMethods/patchDistMethod/patchDistMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchDistMethod, 0);
    defineRunTimeSelectionTable(patchDistMethod, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchDistMethod::patchDistMethod
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    mesh_(mesh),
    patchIDs_(patchIDs)
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::patchDistMethod> Foam::patchDistMethod::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
{
    word patchDistMethodType(dict.lookup("method"));

    Info<< tab << "Selecting patchDistMethod " << patchDistMethodType << endl;

    const auto ctor = ctorTableLookup("patchDistMethodType", dictionaryConstructorTable_(), patchDistMethodType);
    return ctor(dict, mesh, patchIDs);
}



Foam::autoPtr<Foam::patchDistMethod> Foam::patchDistMethod::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const word& defaultPatchDistMethod
)
{
    word patchDistMethodType = defaultPatchDistMethod;
    dict.readIfPresent("method", patchDistMethodType);

    Info<< "Selecting patchDistMethod " << patchDistMethodType << endl;

    const auto ctor = ctorTableLookup("patchDistMethodType", dictionaryConstructorTable_(), patchDistMethodType);
    return ctor(dict, mesh, patchIDs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchDistMethod::~patchDistMethod()
{}


// ************************************************************************* //
