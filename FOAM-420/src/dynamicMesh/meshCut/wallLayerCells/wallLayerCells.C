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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshCut/wallLayerCells/wallLayerCells.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"
#include "meshCut/wallLayerCells/wallNormalInfo/wallNormalInfo.H"
#include "db/IOstreams/Fstreams/OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallLayerCells, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::wallLayerCells::usesCoupledPatch(const label celli) const
{
    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    const cell& cFaces = mesh().cells()[celli];

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        label patchID = patches.whichPatch(facei);

        if ((patchID >= 0) && (patches[patchID].coupled()))
        {
            return true;
        }
    }
    return false;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLayerCells::wallLayerCells
(
    const polyMesh& mesh,
    const List<word>& patchNames,
    const label nLayers
)
:
    edgeVertex(mesh),
    List<refineCell>()
{
    // Find out cells connected to walls.

    const polyPatchList& patches = mesh.boundaryMesh();

    // Make map from name to local patch ID
    HashTable<label> patchNameToIndex(patches.size());

    forAll(patches, patchi)
    {
        patchNameToIndex.insert(patches[patchi].name(), patchi);
    }


    // Count size of walls to set
    label nWalls = 0;

    forAll(patchNames, patchNameI)
    {
        const word& name = patchNames[patchNameI];

        if (patchNameToIndex.found(name))
        {
            label patchi = patchNameToIndex[name];

            nWalls += patches[patchi].size();
        }
    }

    // Allocate storage for start of wave on faces
    List<wallNormalInfo> changedFacesInfo(nWalls);
    labelList changedFaces(nWalls);

    // Fill changedFaces info
    label nChangedFaces = 0;

    forAll(patchNames, patchNameI)
    {
        const word& name = patchNames[patchNameI];

        if (patchNameToIndex.found(name))
        {
            label patchi = patchNameToIndex[name];

            const polyPatch& pp = patches[patchi];

            forAll(pp, patchFacei)
            {
                label meshFacei = pp.start() + patchFacei;

                changedFaces[nChangedFaces] = meshFacei;

                // Set transported information to the wall normal.
                const vector& norm = pp.faceNormals()[patchFacei];

                changedFacesInfo[nChangedFaces] = wallNormalInfo(norm);

                nChangedFaces++;
            }
        }
    }


    // Do a wave of nLayers, transporting the index in patchNames
    // (cannot use local patchIDs since we might get info from neighbouring
    //  processor)

    List<wallNormalInfo> faceInfo(mesh.nFaces()), cellInfo(mesh.nCells());
    FaceCellWave<wallNormalInfo> regionCalc
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceInfo,
        cellInfo,
        0
    );

    regionCalc.iterate(nLayers);


    // Now faceInfo should hold info on faces that are reachable from
    // changedFaces within nLayers iterations. We use face info since that is
    // guaranteed to be consistent across processor boundaries.

    if (debug)
    {
        InfoInFunction << "Dumping selected faces to selectedFaces.obj" << endl;

        OFstream fcStream("selectedFaces.obj");

        label vertI = 0;

        forAll(faceInfo, facei)
        {
            const wallNormalInfo& info = faceInfo[facei];

            if (info.valid(regionCalc.data()))
            {
                const face& f = mesh.faces()[facei];

                point mid(0.0, 0.0, 0.0);

                forAll(f, fp)
                {
                    mid += mesh.points()[f[fp]];
                }
                mid /= f.size();

                fcStream
                    << "v " << mid.x() << ' ' << mid.y() << ' ' << mid.z()
                    << endl;
                vertI++;

                point end(mid + info.normal());

                fcStream
                    << "v " << end.x() << ' ' << end.y() << ' ' << end.z()
                    << endl;
                vertI++;

                fcStream << "l " << vertI << ' ' <<vertI-1 << endl;
            }
        }
    }


    //
    // Copy meshWave information to List<refineCell>
    //

    // Estimate size

    DynamicList<refineCell> refineCells(3*nWalls);

    forAll(cellInfo, celli)
    {
        const wallNormalInfo& info = cellInfo[celli];

        if (info.valid(regionCalc.data()) && !usesCoupledPatch(celli))
        {
            refineCells.append(refineCell(celli, info.normal()));
        }
    }

    // Transfer refineCells storage to this.
    transfer(refineCells);
}


// ************************************************************************* //
