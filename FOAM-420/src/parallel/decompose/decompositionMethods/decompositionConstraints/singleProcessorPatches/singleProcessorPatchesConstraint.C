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
    (c) 2015 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "decompositionConstraints/singleProcessorPatches/singleProcessorPatchesConstraint.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "sets/topoSets/faceSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
defineTypeName(singleProcessorPatchesConstraint);
addToRunTimeSelectionTable
(
    decompositionConstraint,
    singleProcessorPatchesConstraint,
    dictionary
);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::singleProcessorPatchesConstraint::
singleProcessorPatchesConstraint
(
    const dictionary& constraintsDict,
    const word& modelType
)
:
    decompositionConstraint(constraintsDict, typeName),
    patchNamesAndProcs_(coeffDict_.lookup("singleProcessorPatches"))
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep" << endl;

        forAll(patchNamesAndProcs_, setI)
        {
            Info<< "    all cells connected to patchNames "
                << patchNamesAndProcs_[setI].first()
                << " on processor " << patchNamesAndProcs_[setI].second()
                << endl;
        }
    }
}


Foam::decompositionConstraints::singleProcessorPatchesConstraint::
singleProcessorPatchesConstraint
(
    const List<Tuple2<wordReList, label>>& patchNamesAndProcs
)
:
    decompositionConstraint(dictionary(), typeName),
    patchNamesAndProcs_(patchNamesAndProcs)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep" << endl;

        forAll(patchNamesAndProcs_, setI)
        {
            Info<< "    all cells connected to faceSet "
                << patchNamesAndProcs_[setI].first()
                << " on processor " << patchNamesAndProcs_[setI].second()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::singleProcessorPatchesConstraint::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    blockedFace.setSize(mesh.nFaces(), true);

    // Mark faces already in set
    labelList faceToSet(mesh.nFaces(), -1);
    forAll(specifiedProcessorFaces, setI)
    {
        const labelList& faceLabels = specifiedProcessorFaces[setI];
        forAll(faceLabels, i)
        {
            faceToSet[faceLabels[i]] = setI;
        }
    }

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(patchNamesAndProcs_, setI)
    {
        const label destProcI = patchNamesAndProcs_[setI].second();

        const wordReList patchNames = patchNamesAndProcs_[setI].first();
        labelHashSet patchSet = pbm.patchSet(patchNames, false, true);

        labelList nMatch(specifiedProcessorFaces.size(), 0);

        const label nBoundaryFaces = mesh.nFaces() - mesh.nInternalFaces();
        DynamicList<label> markedFaces(nBoundaryFaces);

        forAll(pbm, patchI)
        {
            if (patchSet.found(patchI))
            {
                const polyPatch& pp = pbm[patchI];
                label startI = pp.start();

                forAll(pp, i)
                {
                    label faceI = startI + i;
                    label setI = faceToSet[faceI];
                    markedFaces.append(faceI);

                    if (setI != -1)
                    {
                        nMatch[setI]++;
                    }
                }
            }
        }
        markedFaces.shrink();

        // Only store if all faces are not yet in specifiedProcessorFaces
        // (on all processors)
        bool store = true;

        forAll(nMatch, setI)
        {
            if (nMatch[setI] == markedFaces.size())
            {
                // full match
                store = false;
                break;
            }
            else if (nMatch[setI] > 0)
            {
                // partial match
                store = false;
                break;
            }
        }

        reduce(store, andOp<bool>());

        if (store)
        {
            specifiedProcessorFaces.append(new labelList(markedFaces));
            specifiedProcessor.append(destProcI);
        }
    }

    // Unblock all point connected faces
    // 1. Mark all points on specifiedProcessorFaces
    boolList procFacePoint(mesh.nPoints(), false);
    forAll(specifiedProcessorFaces, setI)
    {
        const labelList& set = specifiedProcessorFaces[setI];
        forAll(set, fI)
        {
            const face& f = mesh.faces()[set[fI]];
            forAll(f, fp)
            {
                procFacePoint[f[fp]] = true;
            }
        }
    }
    syncTools::syncPointList(mesh, procFacePoint, orEqOp<bool>(), false);

    // 2. Unblock all faces on procFacePoint

    label nUnblocked = 0;

    forAll(procFacePoint, pointI)
    {
        if (procFacePoint[pointI])
        {
            const labelList& pFaces = mesh.pointFaces()[pointI];
            forAll(pFaces, i)
            {
                if (blockedFace[pFaces[i]])
                {
                    blockedFace[pFaces[i]] = false;
                    nUnblocked++;
                }
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        reduce(nUnblocked, sumOp<label>());
        Info<< type() << " : unblocked " << nUnblocked << " faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::decompositionConstraints::singleProcessorPatchesConstraint::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // For specifiedProcessorFaces rework the cellToProc to enforce
    // all on one processor since we can't guarantee that the input
    // to regionSplit was a single region.
    // E.g. faceSet 'a' with the cells split into two regions
    // by a notch formed by two walls
    //
    //          \   /
    //           \ /
    //    ---a----+-----a-----
    //
    //
    // Note that reworking the cellToProc might make the decomposition
    // unbalanced.
    label nChanged = 0;

    forAll(specifiedProcessorFaces, setI)
    {
        const labelList& set = specifiedProcessorFaces[setI];

        // Get the processor to use for the set
        label procI = specifiedProcessor[setI];
        if (procI == -1)
        {
            // If no processor specified use the one from the
            // 0th element
            if (set.size())
            {
                procI = decomposition[mesh.faceOwner()[set[0]]];
            }
            reduce(procI, maxOp<label>());
        }

        // Get all points on the sets
        boolList procFacePoint(mesh.nPoints(), false);
        forAll(set, fI)
        {
            const face& f = mesh.faces()[set[fI]];
            forAll(f, fp)
            {
                procFacePoint[f[fp]] = true;
            }
        }
        syncTools::syncPointList(mesh, procFacePoint, orEqOp<bool>(), false);

        // 2. Unblock all faces on procFacePoint
        forAll(procFacePoint, pointI)
        {
            if (procFacePoint[pointI])
            {
                const labelList& pFaces = mesh.pointFaces()[pointI];
                forAll(pFaces, i)
                {
                    label faceI = pFaces[i];

                    label own = mesh.faceOwner()[faceI];
                    if (decomposition[own] != procI)
                    {
                        decomposition[own] = procI;
                        nChanged++;
                    }
                    if (mesh.isInternalFace(faceI))
                    {
                        label nei = mesh.faceNeighbour()[faceI];
                        if (decomposition[nei] != procI)
                        {
                            decomposition[nei] = procI;
                            nChanged++;
                        }
                    }
                }
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        reduce(nChanged, sumOp<label>());
        Info<< type() << " : changed decomposition on " << nChanged
            << " cells" << endl;
    }
}


// ************************************************************************* //
