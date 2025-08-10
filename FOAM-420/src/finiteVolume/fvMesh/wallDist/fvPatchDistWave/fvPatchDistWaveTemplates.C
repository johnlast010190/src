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
void Foam::fvPatchDistWave::setChangedFaces
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    List<labelPair>& changedPatchAndFaces,
    List<PatchPointType>& changedFacesInfo,
    const InitialPatchData& ... initialPatchData
)
{
    label nChangedFaces = 0;
    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        nChangedFaces += mesh.boundary()[iter.key()].size();
    }

    changedPatchAndFaces.resize(nChangedFaces);
    changedFacesInfo.resize(nChangedFaces);

    label changedFacei = 0;

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();

        const fvPatch& patch = mesh.boundary()[patchi];

        forAll(patch.Cf(), patchFacei)
        {
            const label meshFacei = patch.patch().localToGlobal(patchFacei);

            // Filter faces that are counted on direct and indirect patches
            bool countFace = true;
            if (isA<indirectPolyPatch>(patch.patch()))
            {
                // If internal face, add only from master
                if (mesh.isInternalFace(meshFacei))
                {
                    const indirectPolyPatch& ipp =
                        refCast<const indirectPolyPatch>(patch.patch());
                    if (ipp.indirectPolyPatchType() != "master")
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
                changedPatchAndFaces[changedFacei] = {patchi, patchFacei};

                changedFacesInfo[changedFacei] =
                    PatchPointType
                    (
                        patch.Cf()[patchFacei],
                        initialPatchData[patchi][patchFacei] ...,
                        scalar(0)
                    );

                changedFacei++;
            }
        }
    }

    // Resize if needed
    if (changedFacei != changedPatchAndFaces.size())
    {
        changedPatchAndFaces.setSize(changedFacei);
        changedFacesInfo.setSize(changedFacei);
    }
}


template<class PatchPointType, class DataType, class DataMethod>
Foam::label Foam::fvPatchDistWave::getCellValues
(
    FvFaceCellWave<PatchPointType>& waveInfo,
    Field<DataType>& cellValues,
    DataMethod method,
    const DataType& stabiliseValue
)
{
    const List<PatchPointType>& cellInfo = waveInfo.cellInfo();

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


template<class PatchPointType, class DataType, class DataMethod>
Foam::label Foam::fvPatchDistWave::getPatchValues
(
    FvFaceCellWave<PatchPointType>& waveInfo,
    GeometricBoundaryField<DataType, fvPatchField, volMesh>& valuesBf,
    DataMethod method,
    const DataType& stabiliseValue
)
{
    const List<List<PatchPointType>>& patchFaceInfo = waveInfo.patchFaceInfo();

    label nInvalid = 0;

    forAll(valuesBf, patchi)
    {
        forAll(valuesBf[patchi], patchFacei)
        {
            valuesBf[patchi][patchFacei] =
                (patchFaceInfo[patchi][patchFacei].*method)(waveInfo.data())
              + stabiliseValue;

            nInvalid +=
                !patchFaceInfo[patchi][patchFacei].valid(waveInfo.data());
        }
    }

    return nInvalid;
}


template<class PatchPointType>
Foam::label Foam::fvPatchDistWave::wave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    volScalarField& distance,
    const bool correct
)
{
    // Initialise changedFacesInfo to face centres on patches
    List<labelPair> changedPatchAndFaces;
    List<PatchPointType> changedFacesInfo;
    setChangedFaces(mesh, patchIDs, changedPatchAndFaces, changedFacesInfo);

    // Do calculate patch distance by 'growing' from faces
    List<PatchPointType> internalFaceInfo(mesh.nInternalFaces());
    List<List<PatchPointType>> patchFaceInfo
    (
        FvFaceCellWave<PatchPointType>::template
        sizesListList<List<List<PatchPointType>>>
        (
            FvFaceCellWave<PatchPointType>::template
            listListSizes(mesh.boundary()),
            PatchPointType()
        )
    );
    List<PatchPointType> cellInfo(mesh.nCells());

    // Prevent hangs associated with generation of on-demand geometry
    mesh.C();
    mesh.Cf();

    FvFaceCellWave<PatchPointType> wave
    (
        mesh,
        changedPatchAndFaces,
        changedFacesInfo,
        internalFaceInfo,
        patchFaceInfo,
        cellInfo,
        mesh.globalData().nTotalCells() + 1    // max iterations
    );

    // Copy distance into return field
    const label nUnset =
        getCellValues
        (
            wave,
            distance.primitiveFieldRef(),
            &PatchPointType::template dist<int>
        )
      + getPatchValues
        (
            wave,
            distance.boundaryFieldRef(),
            &PatchPointType::template dist<int>,
            SMALL
        );

    // Correct patch cells for true distance
    if (correct)
    {
        Map<labelPair> nearestFace(2*changedFacesInfo.size());
        patchDistFuncs::correctBoundaryFaceCells
        (
            mesh,
            patchIDs,
            distance.primitiveFieldRef(),
            nearestFace
        );
        patchDistFuncs::correctBoundaryPointCells
        (
            mesh,
            patchIDs,
            distance.primitiveFieldRef(),
            nearestFace
        );
    }

    return nUnset;
}


template<class PatchPointType>
Foam::label Foam::fvPatchDistWave::wave
(
    const fvMesh& mesh,
    const labelHashSet& patchIDs,
    const GeometricBoundaryField
        <typename PatchPointType::dataType, fvPatchField, volMesh>&
        initialPatchData,
    volScalarField& distance,
    GeometricField
        <typename PatchPointType::dataType, fvPatchField, volMesh>& data,
    const bool correct
)
{
    // Initialise changedFacesInfo to face centres on patches
    List<labelPair> changedPatchAndFaces;
    List<PatchPointType> changedFacesInfo;
    setChangedFaces
    (
        mesh,
        patchIDs,
        changedPatchAndFaces,
        changedFacesInfo,
        initialPatchData
    );

    // Do calculate patch distance by 'growing' from faces
    List<PatchPointType> internalFaceInfo(mesh.nInternalFaces());
    List<List<PatchPointType>> patchFaceInfo
    (
        FvFaceCellWave<PatchPointType>::template
        sizesListList<List<List<PatchPointType>>>
        (
            FvFaceCellWave<PatchPointType>::template
            listListSizes(mesh.boundary()),
            PatchPointType()
        )
    );
    List<PatchPointType> cellInfo(mesh.nCells());

    // Prevent hangs associated with generation of on-demand geometry
    mesh.C();
    mesh.Cf();

    FvFaceCellWave<PatchPointType> wave
    (
        mesh,
        changedPatchAndFaces,
        changedFacesInfo,
        internalFaceInfo,
        patchFaceInfo,
        cellInfo,
        mesh.globalData().nTotalCells() + 1    // max iterations
    );

    // Copy distance into return field
    const label nUnset =
        getCellValues
        (
            wave,
            distance.primitiveFieldRef(),
            &PatchPointType::template dist<int>
        )
      + getPatchValues
        (
            wave,
            distance.boundaryFieldRef(),
            &PatchPointType::template dist<int>,
            SMALL
        );

    // Copy data into the return field
    getCellValues
    (
        wave,
        data.primitiveFieldRef(),
        &PatchPointType::template data<int>
    );
    getPatchValues
    (
        wave,
        data.boundaryFieldRef(),
        &PatchPointType::template data<int>
    );

    // Correct patch cells for true distance
    if (correct)
    {
        Map<labelPair> nearestPatchAndFace(2*changedFacesInfo.size());
        patchDistFuncs::correctBoundaryFaceCells
        (
            mesh,
            patchIDs,
            distance.primitiveFieldRef(),
            nearestPatchAndFace
        );
        patchDistFuncs::correctBoundaryPointCells
        (
            mesh,
            patchIDs,
            distance.primitiveFieldRef(),
            nearestPatchAndFace
        );

        // Transfer data from nearest face to cell
        forAllConstIter(Map<labelPair>, nearestPatchAndFace, iter)
        {
            const label celli = iter.key();
            const label patchi = iter().first();
            const label patchFacei = iter().second();
            data.primitiveFieldRef()[celli] =
                wave.patchFaceInfo()[patchi][patchFacei].data();
        }
    }

    return nUnset;
}


// ************************************************************************* //
