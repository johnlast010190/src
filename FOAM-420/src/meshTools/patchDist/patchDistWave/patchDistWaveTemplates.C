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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "patchDist/patchDistFuncs/patchDistFuncs.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class PatchPointType, class ... InitialPatchData>
void Foam::patchDistWave::setChangedFaces
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    labelList& changedFaces,
    List<PatchPointType>& changedFacesInfo,
    const InitialPatchData& ... initialPatchData
)
{
    label nChangedFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nChangedFaces += mesh.boundaryMesh()[iter.key()].size();
    }

    changedFaces.resize(nChangedFaces);
    changedFacesInfo.resize(nChangedFaces);

    label changedFacei = 0;

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();

        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        forAll(patch.faceCentres(), patchFacei)
        {
            // const label meshFacei = patch.start() + patchFacei;
            const label meshFacei = patch.localToGlobal(patchFacei);

            // Filter faces that are counted on direct and indirect patches
            bool countFace = true;
            if (isA<indirectPolyPatch>(patch))
            {
                // If internal face, add only from master
                if (mesh.isInternalFace(meshFacei))
                {
                    const indirectPolyPatch& dpp =
                        refCast<const indirectPolyPatch>(patch);
                    if (dpp.indirectPolyPatchType() != "master")
                    {
                        countFace = false;
                    }
                }
                else
                {
                    // Add the face of direct patch only
                    label dPI = mesh.boundaryMesh().whichPatch(meshFacei);
                    if (isA<wall>(mesh.boundaryMesh()[dPI]))
                    {
                        countFace = false;
                    }
                }
            }

            if (countFace)
            {
                changedFaces[changedFacei] = meshFacei;

                changedFacesInfo[changedFacei] =
                    PatchPointType
                    (
                        patch.faceCentres()[patchFacei],
                        initialPatchData[patchi][patchFacei] ...,
                        scalar(0)
                    );

                changedFacei++;
            }
        }
    }

    // Resize if needed
    if (changedFacei != changedFaces.size())
    {
        changedFaces.setSize(changedFacei);
        changedFacesInfo.setSize(changedFacei);
    }
}


template<class PatchPointType, class DataType, class DataMethod>
Foam::label Foam::patchDistWave::getCellValues
(
    FaceCellWave<PatchPointType>& waveInfo,
    Field<DataType>& cellValues,
    DataMethod method,
    const DataType& stabiliseValue
)
{
    const List<PatchPointType>& cellInfo = waveInfo.allCellInfo();

    label nInvalid = 0;

    forAll(cellInfo, celli)
    {
        cellValues[celli] =
            (cellInfo[celli].*method)(waveInfo.data())
          + stabiliseValue;

        nInvalid += !cellInfo[celli].valid(waveInfo.data());
    }

    return nInvalid;
}


template<class PatchPointType>
Foam::label Foam::patchDistWave::wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance,
    const bool correct
)
{
    // Initialise changedFacesInfo to face centres on patches
    List<PatchPointType> changedFacesInfo;
    labelList changedFaces;
    setChangedFaces(mesh, patchIDs, changedFaces, changedFacesInfo);

    // Do calculate patch distance by 'growing' from faces
    List<PatchPointType> faceInfo(mesh.nFaces()), cellInfo(mesh.nCells());
    FaceCellWave<PatchPointType> wave
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceInfo,
        cellInfo,
        mesh.globalData().nTotalCells() + 1    // max iterations
    );

    // Copy distance into return field
    const label nUnset =
        getCellValues
        (
            wave,
            cellDistance,
            &PatchPointType::template dist<int>
        );

    // Correct patch cells for true distance
    if (correct)
    {
        Map<labelPair> nearestFace(2*changedFacesInfo.size());
        patchDistFuncs::correctBoundaryFaceCells
        (
            mesh,
            patchIDs,
            cellDistance,
            nearestFace
        );
        patchDistFuncs::correctBoundaryPointCells
        (
            mesh,
            patchIDs,
            cellDistance,
            nearestFace
        );
    }

    return nUnset;
}


// ************************************************************************* //
