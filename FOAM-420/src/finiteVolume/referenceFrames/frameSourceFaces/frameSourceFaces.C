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
    (c) 2020-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFrames/frameSourceFaces/frameSourceFaces.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "fvMesh/fvPatches/constraint/cyclicACMI/cyclicACMIFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "sets/topoSets/faceSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(frameSourceFaces, 0);
}


// * * * * * * * * * * * *  * * Constructors * * * * *  * * * * * * * * * * //

Foam::frameSourceFaces::frameSourceFaces
(
    const word& name,
    const fvMesh& mesh,
    const labelList& cells,
    const bool& active,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    cells_(cells),
    coeffsDict_(dict),
    nZoneFaces_(0),
    internalFaces_(mesh_.nFaces()),
    includedPatchFaces_(mesh_.boundaryMesh().size()),
    defaultMovingPatchFaces_(mesh_.boundaryMesh().size()),
    faceType_(mesh.nFaces(), notMoving),
    cellInFrame_(mesh_.nCells(), false)
{
    if (active)
    {
        UIndirectList<bool>(cellInFrame_, cells_) = true;
    }
    categorizeFaces();
    createFrameFaceLists();
}


void Foam::frameSourceFaces::categorizeFaces()
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (cellInFrame_[own[facei]] || cellInFrame_[nei[facei]])
        {
            faceType_[facei] = moving;
            nZoneFaces_++;
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (!isA<emptyPolyPatch>(pp))
        {
            if
            (
                pp.coupled()
             || isA<cyclicAMIPolyPatch>(pp)
             || (
                   !includedPatchSet_.empty()
                && !includedPatchSet_.found(patchi)
                )
            )
            {
                forAll(pp, patchFacei)
                {
                    const label facei = pp.start() + patchFacei;
                    if (cellInFrame_[own[facei]])
                    {
                        faceType_[facei] = defaultMoving;
                        nZoneFaces_++;
                    }
                }
            }
            else
            {
                forAll(pp, patchFacei)
                {
                    const label facei = pp.start() + patchFacei;
                    if (cellInFrame_[own[facei]])
                    {
                        faceType_[facei] = moving;
                        nZoneFaces_++;
                    }
                }
            }
        }
    }
    syncTools::syncFaceList(mesh_, faceType_, maxEqOp<label>());
}


void Foam::frameSourceFaces::createFrameFaceLists()
{
    // Set internal faces
    label nInternal = 0;
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceType_[facei] == moving)
        {
            internalFaces_[nInternal++] = facei;
        }
    }
    internalFaces_.setSize(nInternal);


    // Set sub list sizes for included patches
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(patches, patchi)
    {
        label nFacesMoving = 0;
        label nDefaultMovingFaces = 0;
        forAll(patches[patchi], patchFacei)
        {
            const label facei = patches[patchi].start() + patchFacei;
            if (faceType_[facei] == moving)
            {
                nFacesMoving++;
            }
            else if (faceType_[facei] == defaultMoving)
            {
                nDefaultMovingFaces++;
            }
        }
        includedPatchFaces_[patchi].setSize(nFacesMoving);
        defaultMovingPatchFaces_[patchi].setSize(nDefaultMovingFaces);
    }


    // Set face labels for included patches
    forAll(patches, patchi)
    {
        label inFacei = 0;
        label defaultFacei = 0;
        forAll(patches[patchi], patchFacei)
        {
            const label facei = patches[patchi].start() + patchFacei;
            if (faceType_[facei] == moving)
            {
                includedPatchFaces_[patchi][inFacei++] = patchFacei;
            }
            else if (faceType_[facei] == defaultMoving)
            {
                defaultMovingPatchFaces_[patchi][defaultFacei++] = patchFacei;
            }
        }
    }

    if (debug)
    {
        printDebugInformation();
    }
}


void Foam::frameSourceFaces::updateSourceFaces(const labelList& patchIDs)
{
    if (!patchIDs.empty())
    {
        includedPatchSet_ = patchIDs;
        categorizeFaces();
        createFrameFaceLists();
    }
}


const Foam::labelListList& Foam::frameSourceFaces::includedFaces() const
{
    return includedPatchFaces_;
}


const Foam::labelListList& Foam::frameSourceFaces::excludedFaces() const
{
    return defaultMovingPatchFaces_;
}


const Foam::labelList& Foam::frameSourceFaces::internalFaces() const
{
    return internalFaces_;
}


bool Foam::frameSourceFaces::writeData(Ostream& os) const
{
    return true;
}


void Foam::frameSourceFaces::printDebugInformation()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    faceSet internalFaces
    (
        mesh_,
        "internalFaces",
        labelHashSet(internalFaces_)
    );
    Pout<< "Writing " << internalFaces.size()
        << " internal faces in MRF zone to faceSet "
        << internalFaces.name() << endl;
    internalFaces.write();

    faceSet MRFFaces(mesh_, "includedFaces", 100);
    forAll(includedPatchFaces_, patchi)
    {
        forAll(includedPatchFaces_[patchi], i)
        {
            label patchFacei = includedPatchFaces_[patchi][i];
            MRFFaces.insert(patches[patchi].start() + patchFacei);
        }
    }
    Pout<< "Writing " << MRFFaces.size()
        << " patch faces in MRF zone to faceSet "
        << MRFFaces.name() << endl;
    MRFFaces.write();

    faceSet excludedFaces(mesh_, "excludedFaces", 100);
    forAll(defaultMovingPatchFaces_, patchi)
    {
        forAll(defaultMovingPatchFaces_[patchi], i)
        {
            label patchFacei = defaultMovingPatchFaces_[patchi][i];
            excludedFaces.insert(patches[patchi].start() + patchFacei);
        }
    }
    Pout<< "Writing " << excludedFaces.size()
        << " faces in MRF zone with special handling to faceSet "
        << excludedFaces.name() << endl;
    excludedFaces.write();
}


// ************************************************************************* //
