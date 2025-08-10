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

#include "channelIndex.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshTools/meshTools.H"
#include "db/Time/Time.H"
#include "containers/Lists/SortableList/SortableList.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::vector::components,
        3
    >::names[] =
    {
        "x",
        "y",
        "z"
    };
}

const Foam::NamedEnum<Foam::vector::components, 3>
    Foam::channelIndex::vectorComponentsNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Determines face blocking
void Foam::channelIndex::walkOppositeFaces
(
    const polyMesh& mesh,
    const labelList& startFaces,
    boolList& blockedFace
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    DynamicList<label> frontFaces(startFaces);
    forAll(frontFaces, i)
    {
        label facei = frontFaces[i];
        blockedFace[facei] = true;
    }

    while (returnReduce(frontFaces.size(), sumOp<label>()) > 0)
    {
        // Transfer across.
        boolList isFrontBndFace(nBnd, false);
        forAll(frontFaces, i)
        {
            label facei = frontFaces[i];

            if (!mesh.isInternalFace(facei))
            {
                isFrontBndFace[facei-mesh.nInternalFaces()] = true;
            }
        }
        syncTools::swapBoundaryFaceList(mesh, isFrontBndFace);

        // Add
        forAll(isFrontBndFace, i)
        {
            label facei = mesh.nInternalFaces()+i;
            if (isFrontBndFace[i] && !blockedFace[facei])
            {
                blockedFace[facei] = true;
                frontFaces.append(facei);
            }
        }

        // Transfer across cells
        DynamicList<label> newFrontFaces(frontFaces.size());

        forAll(frontFaces, i)
        {
            label facei = frontFaces[i];

            {
                const cell& ownCell = cells[mesh.faceOwner()[facei]];

                label oppositeFacei = ownCell.opposingFaceLabel(facei, faces);

                if (oppositeFacei == -1)
                {
                    FatalErrorInFunction
                        << "Face:" << facei << " owner cell:" << ownCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFacei])
                    {
                        blockedFace[oppositeFacei] = true;
                        newFrontFaces.append(oppositeFacei);
                    }
                }
            }

            if (mesh.isInternalFace(facei))
            {
                const cell& neiCell = mesh.cells()[mesh.faceNeighbour()[facei]];

                label oppositeFacei = neiCell.opposingFaceLabel(facei, faces);

                if (oppositeFacei == -1)
                {
                    FatalErrorInFunction
                        << "Face:" << facei << " neighbour cell:" << neiCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFacei])
                    {
                        blockedFace[oppositeFacei] = true;
                        newFrontFaces.append(oppositeFacei);
                    }
                }
            }
        }

        frontFaces.transfer(newFrontFaces);
    }
}


// Calculate regions.
void Foam::channelIndex::calcLayeredRegions
(
    const polyMesh& mesh,
    const labelList& startFaces
)
{
    boolList blockedFace(mesh.nFaces(), false);
    walkOppositeFaces
    (
        mesh,
        startFaces,
        blockedFace
    );


    if (false)
    {
        OFstream str(mesh.time().path()/"blockedFaces.obj");
        label vertI = 0;
        forAll(blockedFace, facei)
        {
            if (blockedFace[facei])
            {
                const face& f = mesh.faces()[facei];
                forAll(f, fp)
                {
                    meshTools::writeOBJ(str, mesh.points()[f[fp]]);
                }
                str<< 'f';
                forAll(f, fp)
                {
                    str << ' ' << vertI+fp+1;
                }
                str << nl;
                vertI += f.size();
            }
        }
    }


    // Do analysis for connected regions
    cellRegion_.reset(new regionSplit(mesh, blockedFace));

    Info<< "Detected " << cellRegion_().nRegions() << " layers." << nl << endl;

    // Sum number of entries per region
    regionCount_ = regionSum(scalarField(mesh.nCells(), 1.0));

    // Average cell centres to determine ordering.
    pointField regionCc
    (
        regionSum(mesh.cellCentres())
      / regionCount_
    );

    SortableList<scalar> sortComponent(regionCc.component(dir_));

    sortMap_ = sortComponent.indices();

    y_ = sortComponent;

    if (symmetric_)
    {
        y_.setSize(cellRegion_().nRegions()/2);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    symmetric_(readBool(dict.lookup("symmetric"))),
    dir_(vectorComponentsNames_.read(dict.lookup("component")))
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const wordList patchNames(dict.lookup("patches"));

    label nFaces = 0;

    forAll(patchNames, i)
    {
        const label patchi = patches.findPatchID(patchNames[i]);

        if (patchi == -1)
        {
            FatalErrorInFunction
                << "Illegal patch " << patchNames[i]
                << ". Valid patches are " << patches.name()
                << exit(FatalError);
        }

        nFaces += patches[patchi].size();
    }

    labelList startFaces(nFaces);
    nFaces = 0;

    forAll(patchNames, i)
    {
        const polyPatch& pp = patches[patchNames[i]];

        forAll(pp, j)
        {
            startFaces[nFaces++] = pp.start()+j;
        }
    }

    // Calculate regions.
    calcLayeredRegions(mesh, startFaces);
}


Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const labelList& startFaces,
    const bool symmetric,
    const direction dir
)
:
    symmetric_(symmetric),
    dir_(dir)
{
    // Calculate regions.
    calcLayeredRegions(mesh, startFaces);
}


// ************************************************************************* //
