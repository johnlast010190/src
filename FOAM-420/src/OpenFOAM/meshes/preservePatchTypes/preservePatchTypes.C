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
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/preservePatchTypes/preservePatchTypes.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMeshEntries.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::preservePatchTypes
(
    const objectRegistry& obr,
    const word& meshInstance,
    const fileName& meshDir,
    const wordList& patchNames,
    PtrList<dictionary>& patchDicts,
    const word& defaultFacesName,
    word& defaultFacesType
)
{
    patchDicts.setSize(patchNames.size());

    dictionary patchDictionary;

    // Read boundary file as single dictionary
    {
        IOobject patchEntriesHeader
        (
            "boundary",
            meshInstance,
            meshDir,
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (patchEntriesHeader.typeHeaderOk<polyBoundaryMesh>(true))
        {
            // Create a list of entries from the boundary file.
            polyBoundaryMeshEntries patchEntries(patchEntriesHeader);

            forAll(patchEntries, patchi)
            {
                patchDictionary.add(patchEntries[patchi]);
            }
        }
    }

    forAll(patchNames, patchi)
    {
        if (patchDictionary.found(patchNames[patchi]))
        {
            const dictionary& patchDict =
                patchDictionary.subDict(patchNames[patchi]);

            patchDicts.set(patchi, patchDict.clone());
            patchDicts[patchi].remove("nFaces");
            patchDicts[patchi].remove("startFace");
        }
    }

    if (patchDictionary.found(defaultFacesName))
    {
        const dictionary& patchDict =
            patchDictionary.subDict(defaultFacesName);

        patchDict.readIfPresent("geometricType", defaultFacesType);
    }

    Info<< nl << "Default patch type set to " << defaultFacesType << endl;
}


// ************************************************************************* //
