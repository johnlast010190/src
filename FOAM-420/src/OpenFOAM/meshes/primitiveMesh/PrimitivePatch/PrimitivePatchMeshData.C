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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
#include "containers/HashTables/Map/Map.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcMeshData() const
{
    if (debug)
    {
        PoutInFunction << "Calculating mesh data in PrimitivePatch" << endl;
    }

    // It is considered an error to attempt to recalculate meshPoints
    // if they have already been calculated.
    if (meshPointsPtr_ || localFacesPtr_)
    {
        FatalErrorInFunction
            << "meshPointsPtr_ or localFacesPtr_ already allocated"
            << abort(FatalError);
    }

    // Create a map for marking points.  Estimated size is 4 times the
    // number of faces in the patch
    Map<label> markedPoints(4*this->size());


    // Important:
    // ~~~~~~~~~~
    // In <= 1.5 the meshPoints would be in increasing order but this gives
    // problems in processor point synchronisation where we have to find out
    // how the opposite side would have allocated points.

    ////- 1.5 code:
    //// if the point is used, set the mark to 1
    //forAll(*this, facei)
    //{
    //    const FaceType& curPoints = this->operator[](facei);
    //
    //    forAll(curPoints, pointi)
    //    {
    //        markedPoints.insert(curPoints[pointi], -1);
    //    }
    //}
    //
    //// Create the storage and store the meshPoints.  Mesh points are
    //// the ones marked by the usage loop above
    //meshPointsPtr_ = new labelList(markedPoints.toc());
    //labelList& pointPatch = *meshPointsPtr_;
    //
    //// Sort the list to preserve compatibility with the old ordering
    //sort(pointPatch);
    //
    //// For every point in map give it its label in mesh points
    //forAll(pointPatch, pointi)
    //{
    //    markedPoints.find(pointPatch[pointi])() = pointi;
    //}

    //- Unsorted version:
    DynamicList<label> meshPoints(2*this->size());
    forAll(*this, facei)
    {
        const FaceType& curPoints = this->operator[](facei);

        forAll(curPoints, pointi)
        {
            if (markedPoints.insert(curPoints[pointi], meshPoints.size()))
            {
                meshPoints.append(curPoints[pointi]);
            }
        }
    }
    // Transfer to straight list (reuses storage)
    meshPointsPtr_ = new labelList(meshPoints, true);


    // Create local faces. Note that we start off from copy of original face
    // list (even though vertices are overwritten below). This is done so
    // additional data gets copied (e.g. region number of labelledTri)
    localFacesPtr_ = new List<FaceType>(*this);
    List<FaceType>& lf = *localFacesPtr_;

    forAll(*this, facei)
    {
        const FaceType& curFace = this->operator[](facei);
        lf[facei].setSize(curFace.size());

        forAll(curFace, labelI)
        {
            lf[facei][labelI] = markedPoints.find(curFace[labelI])();
        }
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating mesh data in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcMeshPointMap() const
{
    if (debug)
    {
        PoutInFunction << "Calculating mesh point map in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate meshPoints
    // if they have already been calculated.
    if (meshPointMapPtr_)
    {
        FatalErrorInFunction
            << "meshPointMapPtr_ already allocated"
            << abort(FatalError);
    }

    const labelList& mp = meshPoints();

    meshPointMapPtr_ = new Map<label>(2*mp.size());
    Map<label>& mpMap = *meshPointMapPtr_;

    forAll(mp, i)
    {
        mpMap.insert(mp[i], i);
    }

    if (debug)
    {
        PoutInFunction
            << "Finished calculating mesh point map in PrimitivePatch" << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcLocalPoints() const
{
    if (debug)
    {
        PoutInFunction << "Calculating localPoints in PrimitivePatch" << endl;
    }

    // It is considered an error to attempt to recalculate localPoints
    // if they have already been calculated.
    if (localPointsPtr_)
    {
        FatalErrorInFunction
            << "localPointsPtr_ already allocated"
            << abort(FatalError);
    }

    const labelList& meshPts = meshPoints();

    localPointsPtr_ = new Field<PointType>(meshPts.size());

    Field<PointType>& locPts = *localPointsPtr_;

    forAll(meshPts, pointi)
    {
        locPts[pointi] = points_[meshPts[pointi]];
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating localPoints in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcPointNormals() const
{
    if (debug)
    {
        PoutInFunction << "Calculating pointNormals in PrimitivePatch" << endl;
    }

    // It is considered an error to attempt to recalculate pointNormals
    // if they have already been calculated.
    if (pointNormalsPtr_)
    {
        FatalErrorInFunction
            << "pointNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    const Field<PointType>& faceUnitNormals = faceNormals();

    const labelListList& pf = pointFaces();

    pointNormalsPtr_ = new Field<PointType>
    (
        meshPoints().size(),
        PointType::zero
    );

    Field<PointType>& n = *pointNormalsPtr_;

    forAll(pf, pointi)
    {
        PointType& curNormal = n[pointi];

        const labelList& curFaces = pf[pointi];

        forAll(curFaces, facei)
        {
            curNormal += faceUnitNormals[curFaces[facei]];
        }

    if (mag(curNormal) < VSMALL)
    {
        if (curFaces.size())
        {
            curNormal = faceUnitNormals[curFaces[0]];
        }
        else
        {
            curNormal = vector(1., 0, 0);
        }
    }
    else
    {
        curNormal /= mag(curNormal) + VSMALL;
    }
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating pointNormals in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcFaceCentres() const
{
    if (debug)
    {
        PoutInFunction << "Calculating faceCentres in PrimitivePatch" << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (faceCentresPtr_)
    {
        FatalErrorInFunction
            << "faceCentresPtr_ already allocated"
            << abort(FatalError);
    }

    faceCentresPtr_ = new Field<PointType>(this->size());

    Field<PointType>& c = *faceCentresPtr_;

    forAll(c, facei)
    {
        c[facei] = this->operator[](facei).centre(points_);
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating faceCentres in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcMagFaceAreas() const
{
    if (debug)
    {
        PoutInFunction << "Calculating magFaceAreas in PrimitivePatch" << endl;
    }

    // It is an error to calculate these more than once
    if (magFaceAreasPtr_)
    {
        FatalErrorInFunction
            << "magFaceAreasPtr_ already allocated"
            << abort(FatalError);
    }

    magFaceAreasPtr_ = new Field<scalar>(this->size());
    Field<scalar>& a = *magFaceAreasPtr_;

    forAll(a, facei)
    {
        a[facei] = this->operator[](facei).mag(points_);
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating magFaceAreas in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcFaceAreas() const
{
    if (debug)
    {
        PoutInFunction << "Calculating faceAreas in PrimitivePatch" << endl;
    }

    // It is considered an error to attempt to recalculate faceNormals
    // if they have already been calculated.
    if (faceAreasPtr_)
    {
        FatalErrorInFunction
            << "faceAreasPtr_ already allocated"
            << abort(FatalError);
    }

    faceAreasPtr_ = new Field<PointType>(this->size());

    Field<PointType>& n = *faceAreasPtr_;

    forAll(n, facei)
    {
        n[facei] = this->operator[](facei).areaNormal(points_);
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating faceAreas in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcFaceNormals() const
{
    if (debug)
    {
        PoutInFunction << "Calculating faceNormals in PrimitivePatch" << endl;
    }

    // It is considered an error to attempt to recalculate faceNormals
    // if they have already been calculated.
    if (faceNormalsPtr_)
    {
        FatalErrorInFunction
            << "faceNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    faceNormalsPtr_ = new Field<PointType>(this->size());

    Field<PointType>& n = *faceNormalsPtr_;

    forAll(n, facei)
    {
        n[facei] = this->operator[](facei).unitNormal(points_);
    }

    if (debug)
    {
        PoutInFunction << "Finished calculating faceNormals in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::correctZeroSizedFaceNormals() const
{
    if (debug)
    {
        PoutInFunction << "Correcting faceNormals in PrimitivePatch" << endl;
    }

    Field<PointType>& fn = *faceNormalsPtr_;
    const Field<PointType>& pn = pointNormals();

    forAll(fn, facei)
    {
        if (mag(fn[facei]) == 0)
        {
            //const labelList& curFaces = pf[pointi];
            const FaceType& curFace = localFaces()[facei];
            forAll(curFace, pI)
            {
                fn[facei] += pn[curFace[pI]];
            }
        }
    }

    if (debug)
    {
        PoutInFunction << "Finished correcting faceNormals in PrimitivePatch"
            << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcFaceCurvature() const
{
    if (debug)
    {
        PoutInFunction << "Calculating faceCurvature in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCurvature
    // if they have already been calculated.
    if (faceCurvaturePtr_)
    {
        FatalErrorInFunction
            << "faceCurvaturePtr_already allocated"
            << abort(FatalError);
    }

    faceCurvaturePtr_ = new Field<scalar>(this->size());

    Field<scalar>& c = *faceCurvaturePtr_;
    const Field<PointType>& n = pointNormals();

    forAll(c, faceI)
    {
        c[faceI] = 0.;
        PointType nf = this->operator[](faceI).unitNormal(points_);
        PointType fc = this->operator[](faceI).centre(points_ );

        const FaceType& curFace = localFaces()[faceI];

        forAll(curFace, labelI)
        {
            label pointI = curFace[labelI];
            PointType nI = n[pointI];
            PointType curvDir = nI - nf;
            scalar magDir = mag(curvDir);

            if (magDir > SMALL)
            {
                curvDir /=  magDir + VSMALL;

        scalar minDx = 0.1*mag(points_[pointI] - fc);
                scalar dx = fabs((points_[pointI] - fc) & curvDir);
        dx = max(dx, minDx);
                if (dx > SMALL)
                {
                    scalar norm = nI & nf;
                    norm = acos(min(max(norm,-1.),1.));
                    norm /= dx;
                    c[faceI] = max(c[faceI], norm);
                }
            }
        }
    }

    if (debug)
    {
        PoutInFunction
            << "Finished calculating faceCurvature in PrimitivePatch" << endl;
    }
}


template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcPointCurvature
(
    const boolList& flipMap
) const
{
    if (debug)
    {
        PoutInFunction << "Calculating faceCurvature in PrimitivePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCurvature
    // if they have already been calculated.
    if (pointCurvaturePtr_)
    {
        FatalErrorInFunction
            << "faceCurvaturePtr_already allocated"
            << abort(FatalError);
    }

    pointCurvaturePtr_ = new Field<scalar>(localPoints().size());
    Field<scalar>& pCurv = *pointCurvaturePtr_;

    const Field<PointType>& pts = localPoints();
    const labelListList& pFaces = pointFaces();
    const Field<PointType>& fNormals =  faceNormals();
    const Field<PointType>& cf = faceCentres();
    Field<PointType> pointNormal =
        Field<PointType>(pts.size(), PointType::zero);

    forAll(pts, pI)
    {
        const labelList& pFacesI = pFaces[pI];
        PointType normalSum = PointType::zero;
        forAll(pFacesI, fI)
        {
            label fII = pFacesI[fI];
            if (!flipMap[fII])
            {
                normalSum += fNormals[fII];
            }
            else
            {
                normalSum -= fNormals[fII];
            }
        }
        pointNormal[pI] = normalSum/pFacesI.size();
    }

    forAll(pts, pI)
    {
        const labelList& pFacesI = pFaces[pI];
        scalar normalSum = scalar(0.0);
        forAll(pFacesI, fI)
        {
            label fII = pFacesI[fI];

            PointType difN = pointNormal[pI] - fNormals[fII];
            PointType difCoord = pts[pI] - cf[fII];

            scalar lCurv = difN&difCoord;
            normalSum += lCurv;
        }
        pCurv[pI] = normalSum/pFacesI.size();
    }

    if (debug)
    {
        PoutInFunction
            << "Finished calculating faceCurvature in PrimitivePatch" << endl;
    }
}


// ************************************************************************* //
