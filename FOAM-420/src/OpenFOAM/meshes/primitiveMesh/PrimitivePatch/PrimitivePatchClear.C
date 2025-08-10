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
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
#include "include/demandDrivenData.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::clearGeom()
{
    if (debug)
    {
        InfoInFunction << "Clearing geometric data" << endl;
    }

    deleteDemandDrivenData(localPointsPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
    deleteDemandDrivenData(magFaceAreasPtr_);
    deleteDemandDrivenData(faceNormalsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
    deleteDemandDrivenData(faceCurvaturePtr_);
    deleteDemandDrivenData(pointCurvaturePtr_);
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::clearTopology()
{
    if (debug)
    {
        InfoInFunction << "Clearing patch addressing" << endl;
    }

    // group created and destroyed together
    if (edgesPtr_ && faceFacesPtr_ && edgeFacesPtr_ && faceEdgesPtr_)
    {
        delete edgesPtr_;
        edgesPtr_ = nullptr;

        delete faceFacesPtr_;
        faceFacesPtr_ = nullptr;

        delete edgeFacesPtr_;
        edgeFacesPtr_ = nullptr;

        delete faceEdgesPtr_;
        faceEdgesPtr_ = nullptr;
    }

    deleteDemandDrivenData(boundaryPointsPtr_);
    deleteDemandDrivenData(pointEdgesPtr_);
    deleteDemandDrivenData(pointFacesPtr_);
    deleteDemandDrivenData(edgeLoopsPtr_);
    deleteDemandDrivenData(localPointOrderPtr_);
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::clearPatchMeshAddr()
{
    if (debug)
    {
        InfoInFunction << "Clearing patch-mesh addressing" << endl;
    }

    deleteDemandDrivenData(meshPointsPtr_);
    deleteDemandDrivenData(meshPointMapPtr_);
    deleteDemandDrivenData(localFacesPtr_);
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::clearMeshPoints()
{
    if (meshPointsPtr_)
    {
        delete meshPointsPtr_;
        meshPointsPtr_ = nullptr;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::clearOut()
{
    clearGeom();
    clearTopology();
    clearPatchMeshAddr();
    clearMeshPoints();
}


// ************************************************************************* //
