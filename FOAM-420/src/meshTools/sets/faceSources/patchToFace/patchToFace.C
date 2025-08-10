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

#include "sets/faceSources/patchToFace/patchToFace.H"
#include "meshes/polyMesh/polyMesh.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(patchToFace, 0);

addToRunTimeSelectionTable(topoSetSource, patchToFace, word);

addToRunTimeSelectionTable(topoSetSource, patchToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::patchToFace::usage_
(
    patchToFace::typeName,
    "\n    Usage: patchToFace patch\n\n"
    "    Select all faces in the patch. Note:accepts wildcards for patch.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchToFace::combine(topoSet& set, const bool add) const
{
    labelHashSet patchIDs = mesh_.boundaryMesh().patchSet
    (
        List<wordRe>(1, patchName_),
        true,           // warn if not found
        true            // use patch groups if available
    );

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        label patchi = iter.key();

        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        Info<< "    Found matching patch " << pp.name()
            << " with " << returnReduce(pp.size(), sumOp<label>())
            << " faces." << endl;

        for
        (
            label facei = pp.start();
            facei < pp.start() + pp.size();
            facei++
        )
        {
            addOrDelete(set, facei, add);
        }
    }

    if (patchIDs.empty())
    {
        WarningInFunction
            << "Cannot find any patch named " << patchName_ << endl
            << "Valid names are " << mesh_.boundaryMesh().names() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::patchToFace::patchToFace
(
    const polyMesh& mesh,
    const word& patchName
)
:
    topoSetSource(mesh),
    patchName_(patchName)
{}


// Construct from dictionary
Foam::patchToFace::patchToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    patchName_(dict.lookup("name"))
{}


// Construct from Istream
Foam::patchToFace::patchToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    patchName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToFace::~patchToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all faces of patch " << patchName_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all faces of patch " << patchName_ << " ..."
            << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
