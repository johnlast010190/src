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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2008 OpenCFD Ltd.

Description

\*---------------------------------------------------------------------------*/

#include "sets/cellSources/patchToCell/patchToCell.H"
#include "meshes/polyMesh/polyMesh.H"
#include "sets/topoSets/faceSet.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(patchToCell, 0);

addToRunTimeSelectionTable(topoSetSource, patchToCell, word);

addToRunTimeSelectionTable(topoSetSource, patchToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::patchToCell::usage_
(
    patchToCell::typeName,
    "\n    Usage: patchToCell (<patch0> <patch1> .. <patchn>)\n\n"
    "    Select cells that are next to the given pacthes\n\n"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToCell::combine(topoSet& set, const bool add) const
{
    forAll(patchNames_, pnI)
    {
        label patchI = mesh_.boundaryMesh().findPatchID(patchNames_[pnI]);

        if (patchI != -1)
        {
            forAll(mesh_.boundaryMesh()[patchI], pfI)
            {
                label faceI = mesh_.boundaryMesh()[patchI].start() + pfI;
                label cellI = mesh_.faceOwner()[faceI];
                addOrDelete(set, cellI, add);
            }
        }
        else
        {
            WarningInFunction
                << "Specified patch name not found for cellSet source."
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::patchToCell::patchToCell
(
    const polyMesh& mesh,
    const wordList& patchNames
)
:
    topoSetSource(mesh),
    patchNames_(patchNames)
{}


// Construct from dictionary
Foam::patchToCell::patchToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    patchNames_(dict.lookup("patches"))
{}


// Construct from Istream
Foam::patchToCell::patchToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    patchNames_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToCell::~patchToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells adjacent to patch(es) " << patchNames_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells according to patch(es) " << patchNames_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
