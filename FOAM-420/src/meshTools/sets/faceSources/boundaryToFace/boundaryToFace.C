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

#include "sets/faceSources/boundaryToFace/boundaryToFace.H"
#include "meshes/polyMesh/polyMesh.H"

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boundaryToFace, 0);

addToRunTimeSelectionTable(topoSetSource, boundaryToFace, word);

addToRunTimeSelectionTable(topoSetSource, boundaryToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::boundaryToFace::usage_
(
    boundaryToFace::typeName,
    "\n    Usage: boundaryToFace\n\n"
    "    Select all boundary faces\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundaryToFace::combine(topoSet& set, const bool add) const
{
    for
    (
        label facei = mesh().nInternalFaces();
        facei < mesh().nFaces();
        facei++
    )
    {
        addOrDelete(set, facei, add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::boundaryToFace::boundaryToFace(const polyMesh& mesh)
:
    topoSetSource(mesh)
{}


// Construct from dictionary
Foam::boundaryToFace::boundaryToFace(const polyMesh& mesh, const dictionary&)
:
    topoSetSource(mesh)
{}


// Construct from Istream
Foam::boundaryToFace::boundaryToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryToFace::~boundaryToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all boundary faces ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all boundary faces ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
