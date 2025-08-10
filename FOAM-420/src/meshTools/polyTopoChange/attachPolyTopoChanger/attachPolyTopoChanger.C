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

#include "polyTopoChange/attachPolyTopoChanger/attachPolyTopoChanger.H"
#include "meshes/polyMesh/polyMesh.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::attachPolyTopoChanger::attachPolyTopoChanger
(
    const IOobject& io,
    polyMesh& mesh
)
:
    polyTopoChanger(io, mesh)
{}


Foam::attachPolyTopoChanger::attachPolyTopoChanger
(
    polyMesh& mesh
)
:
    polyTopoChanger(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::attachPolyTopoChanger::attach(const bool removeEmptyPatches)
{
    if (debug)
    {
        Pout<< "void attachPolyTopoChanger::attach(): "
            << "Attaching mesh" << endl;
    }

    // Save current file instance
    const fileName oldInst = mesh_.facesInstance();

    // Execute all polyMeshModifiers
    changeMesh(false);  // no inflation

    const pointField p = mesh_.oldPoints();

    mesh_.movePoints(p);

    if (debug)
    {
        Pout<< "Clearing mesh." << endl;
    }

    if (removeEmptyPatches)
    {
        // Re-do the boundary patches, removing the ones with zero size
        const polyBoundaryMesh& oldPatches = mesh_.boundaryMesh();

        List<polyPatch*> newPatches(oldPatches.size());
        label nNewPatches = 0;

        forAll(oldPatches, patchi)
        {
            if (oldPatches[patchi].size())
            {
                const directPolyPatch& dpp =
                    refCast<const directPolyPatch>(oldPatches[patchi]);
                newPatches[nNewPatches] = dpp.clone
                (
                    mesh_.boundaryMesh(),
                    nNewPatches,
                    oldPatches[patchi].size(),
                    oldPatches[patchi].start()
                ).ptr();

                nNewPatches++;
            }
            else
            {
                if (debug)
                {
                    Pout<< "Removing zero-sized patch " << patchi
                        << " named " << oldPatches[patchi].name() << endl;
                }
            }
        }

        newPatches.setSize(nNewPatches);

        mesh_.removeBoundary();
        mesh_.addPatches(newPatches);
    }

    // Reset the file instance to overwrite the original mesh
    mesh_.setInstance(oldInst);

    if (debug)
    {
        Pout<< "void attachPolyTopoChanger::attach(): "
            << "Finished attaching mesh" << endl;
    }

    mesh_.checkMesh();
}


// ************************************************************************* //
