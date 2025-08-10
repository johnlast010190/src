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
    (c) 1991-2008 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

Description
    Writes out simple data to VTK format for visualisation

\*----------------------------------------------------------------------------*/

#include "meshes/meshTools/simpleVTKWriter.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/boundBox/boundBox.H"

Foam::simpleVTKWriter::simpleVTKWriter() :
    canAddPointData_(true),
    canAddFaceData_(true)
{}

Foam::simpleVTKWriter::simpleVTKWriter
(
    const List<face> &faces,
    const pointField &points,
    const bool removeUnusedPoints
)
:
    canAddPointData_(true),
    canAddFaceData_(true)
{
    if (removeUnusedPoints)
    {
        // Process the data so that unused points aren't written to file
        pointMap_.setSize(points.size(), -1);
        faces_.setSize(faces.size());

        label ptCount = 0;
        forAll(faces, fI)
        {
            faces_[fI].setSize(faces[fI].size());
            forAll(faces[fI],pI)
            {
                if (pointMap_[faces[fI][pI]] == -1)
                {
                    pointMap_[faces[fI][pI]] = ptCount++;
                }
                faces_[fI][pI] = pointMap_[faces[fI][pI]];
            }
        }

        points_.setSize(ptCount);

        forAll(pointMap_,pI)
        {
            if (pointMap_[pI] != -1)
            {
                points_[pointMap_[pI]] = points[pI];
            }
        }
    }
    else
    {
        points_ = points;
        faces_ = faces;
    }
}

Foam::simpleVTKWriter::simpleVTKWriter
(
    const List<face> &faces,
    const labelList &faceList,
    const pointField &points
)
:
    canAddPointData_(true),
    canAddFaceData_(true)
{
    //process the data so that unused points aren't written to file

    pointMap_.setSize(points.size(),-1);
    faceMap_.setSize(faces.size(),-1);

    faces_.setSize(faceList.size());

    label ptCount = 0;
    forAll(faceList,i)
    {
        label fI = faceList[i];
        faceMap_[fI] = i;
        faces_[i].setSize(faces[fI].size());
        forAll(faces[fI],pI)
        {
            if (pointMap_[faces[fI][pI]] == -1)
            {
                pointMap_[faces[fI][pI]] = ptCount++;
            }
            faces_[i][pI] = pointMap_[faces[fI][pI]];
        }
    }

    points_.setSize(ptCount);

    forAll(pointMap_,pI)
    {
        if (pointMap_[pI] != -1)
        {
            points_[pointMap_[pI]] = points[pI];
        }
    }
}


// Really simple version for writing out a line defined as a list points
Foam::simpleVTKWriter::simpleVTKWriter
(
    const labelList &pointList,
    const List<vector> &points
)
:
    canAddPointData_(true),
    canAddFaceData_(false)
{
    if (pointList.size() > 0)
    {
        points_.setSize(pointList.size());
        forAll(points_,j) {
            points_[j] = points[pointList[j]];
        }

        faces_.setSize(points_.size()-1);
        forAll(faces_,j)
        {
            faces_[j].setSize(3);
            faces_[j][0] = j;
            faces_[j][1] = j+1;
            faces_[j][2] = j;
        }
    }
}


// Really simple version for writing out a set of lines defined as lists of points. Doesn't compact numbering.
Foam::simpleVTKWriter::simpleVTKWriter
(
    const List<labelList> &pointList,
    const List<vector> &points
)
:
    points_(points),
    canAddPointData_(true),
    canAddFaceData_(false)
{
    label numFaces = 0;
    forAll(pointList,i)
    {
        numFaces += pointList[i].size()-1;
    }

    if (numFaces > 0)
    {
        faces_.setSize(numFaces);
        label fI = 0;
        forAll(pointList,i)
        {
            for (label j =0; j < pointList[i].size()-1; ++j) {
                faces_[fI].setSize(3);
                faces_[fI][0] = pointList[i][j];
                faces_[fI][1] = pointList[i][j+1];
                faces_[fI][2] = pointList[i][j];
                fI++;
            }
        }
    }
}


// Really simple version for writing out a line defined as a list points
Foam::simpleVTKWriter::simpleVTKWriter(const List<vector> &points)
:
    points_(points),
    canAddPointData_(true),
    canAddFaceData_(false)
{
    if (points.size() > 0)
    {
        faces_.setSize(points_.size()-1);
        forAll(faces_,j)
        {
            faces_[j].setSize(3);
            faces_[j][0] = j;
            faces_[j][1] = j+1;
            faces_[j][2] = j;
        }
    }
}

Foam::simpleVTKWriter::simpleVTKWriter
(
    const List<edge> &edges,
    const List<vector> &points
)
:
    points_(points),
    canAddPointData_(true),
    canAddFaceData_(false)
{
    pointMap_ = identity(points.size());
    faces_.setSize(edges.size());
    forAll(edges,j)
    {
        faces_[j].setSize(3);
        faces_[j][0] = edges[j][0];
        faces_[j][1] = edges[j][1];
        faces_[j][2] = edges[j][0];
    }
}

// Really simple version for writing out a line defined as a list points
Foam::simpleVTKWriter::simpleVTKWriter
(
    const List<edge> &edges,
    const List<vector> &points,
    const DynamicList<label> &edgeList
)
:
    points_(points),
    canAddPointData_(true),
    canAddFaceData_(false)
{
    faces_.setSize(edgeList.size());
    forAll(edgeList,j)
    {
        faces_[j].setSize(3);
        faces_[j][0] = edges[edgeList[j]][0];
        faces_[j][1] = edges[edgeList[j]][1];
        faces_[j][2] = edges[edgeList[j]][0];
    }
}

Foam::simpleVTKWriter::simpleVTKWriter(const boundBox &bbox) :
    canAddPointData_(true),
    canAddFaceData_(true)
{
    points_.setSize(8);

    points_[0] = bbox.min();
    points_[1] = vector(bbox.max().x(),bbox.min().y(),bbox.min().z());
    points_[2] = vector(bbox.max().x(),bbox.max().y(),bbox.min().z());
    points_[3] = vector(bbox.min().x(),bbox.max().y(),bbox.min().z());
    points_[4] = vector(bbox.min().x(),bbox.min().y(),bbox.max().z());
    points_[5] = vector(bbox.max().x(),bbox.min().y(),bbox.max().z());
    points_[6] = bbox.max();
    points_[7] = vector(bbox.min().x(),bbox.max().y(),bbox.max().z());

    faces_.setSize(6);

    faces_[0].setSize(4);
    faces_[0][0] = 0;
    faces_[0][1] = 1;
    faces_[0][2] = 2;
    faces_[0][3] = 3;

    faces_[1].setSize(4);
    faces_[1][0] = 0;
    faces_[1][1] = 1;
    faces_[1][2] = 5;
    faces_[1][3] = 4;

    faces_[2].setSize(4);
    faces_[2][0] = 1;
    faces_[2][1] = 2;
    faces_[2][2] = 6;
    faces_[2][3] = 5;

    faces_[3].setSize(4);
    faces_[3][0] = 2;
    faces_[3][1] = 3;
    faces_[3][2] = 7;
    faces_[3][3] = 6;

    faces_[4].setSize(4);
    faces_[4][0] = 3;
    faces_[4][1] = 0;
    faces_[4][2] = 4;
    faces_[4][3] = 7;

    faces_[5].setSize(4);
    faces_[5][0] = 4;
    faces_[5][1] = 5;
    faces_[5][2] = 6;
    faces_[5][3] = 7;
}


void Foam::simpleVTKWriter::addMorePoints(const List<vector> &points)
{
    if (points.size() > 0)
    {
        //append new points
        label poffset = points_.size();
        points_.setSize(points_.size() + points.size());
        for (label j = poffset; j < points_.size(); ++j)
        {
            points_[j] = points[j-poffset];
        }

        //append new "faces"
        label foffset = faces_.size();
        label diffoffset = poffset - foffset;
        faces_.setSize(faces_.size() + points.size()-1);
        for (label j = foffset; j < faces_.size(); ++j)
        {
            faces_[j].setSize(3);
            faces_[j][0] = j + diffoffset;
            faces_[j][1] = j+1 + diffoffset;
            faces_[j][2] = j + diffoffset;
        }
    }
}

void Foam::simpleVTKWriter::addMorePoints
(
    const labelList &pointList,
    const List<vector> &points
)
{
    if (pointList.size() > 0)
    {
        //append new points
        label poffset = points_.size();
        points_.setSize(points_.size() + pointList.size());
        for (label j = poffset; j < points_.size(); ++j)
        {
            points_[j] = points[pointList[j-poffset]];
        }

        //append new "faces"
        label foffset = faces_.size();
        label diffoffset = poffset - foffset;
        faces_.setSize(faces_.size() + pointList.size()-1);
        for (label j = foffset; j < faces_.size(); ++j)
        {
            faces_[j].setSize(3);
            faces_[j][0] = j + diffoffset;
            faces_[j][1] = j+1 + diffoffset;
            faces_[j][2] = j + diffoffset;
        }
    }
}


void Foam::simpleVTKWriter::addFaceData
(
    const word &name,
    const List<label> &faceData
)
{
    if (canAddFaceData_)
    {
        faceLabelNames_.append(name);
        if (faceMap_.empty())
        {
            faceLabelData_.append(faceData);
        }
        else
        {
            List<label> newFaceData(faces_.size());
            forAll(faceMap_,i)
            {
                if (faceMap_[i] != -1)
                {
                    newFaceData[faceMap_[i]] = faceData[i];
                }
            }
            faceLabelData_.append(newFaceData);
        }
    }
}



void Foam::simpleVTKWriter::addFaceData
(
    const word &name,
    const List<scalar> &faceData
)
{
    if (canAddFaceData_)
    {
        faceScalarNames_.append(name);
        if (faceMap_.empty())
        {
            faceScalarData_.append(faceData);
        }
        else
        {
            List<scalar> newFaceData(faces_.size());
            forAll(faceMap_,i)
            {
                if (faceMap_[i] != -1)
                {
                    newFaceData[faceMap_[i]] = faceData[i];
                }
            }
            faceScalarData_.append(newFaceData);
        }
    }
}

void Foam::simpleVTKWriter::addFaceData
(
    const word &name,
    const List<vector> &faceData
)
{
    if (canAddFaceData_)
    {
        faceVectorNames_.append(name);
        if (faceMap_.empty())
        {
            faceVectorData_.append(faceData);
        }
        else
        {
            List<vector> newFaceData(faces_.size());
            forAll(faceMap_,i)
            {
                if (faceMap_[i] != -1)
                {
                    newFaceData[faceMap_[i]] = faceData[i];
                }
            }
            faceVectorData_.append(newFaceData);
        }
    }
}


void Foam::simpleVTKWriter::addPointData
(
    const word &name,
    const List<label> &pointData
)
{
    if (canAddPointData_)
    {
        pointLabelNames_.append(name);
        pointLabelData_.append(List<label>(points_.size()));
        List<label> &data = pointLabelData_[pointLabelData_.size()-1];
        forAll(pointData,pI)
        {
            if (pointMap_[pI] != -1)
            {
                data[pointMap_[pI]] = pointData[pI];
            }
        }
    }
}


void Foam::simpleVTKWriter::addPointData
(
    const word &name,
    const List<scalar> &pointData
)
{
    if (canAddPointData_)
    {
        pointScalarNames_.append(name);
        pointScalarData_.append(List<scalar>(points_.size()));

        List<scalar> &data = pointScalarData_[pointScalarData_.size()-1];

        forAll(pointData,pI)
        {
            if (pointMap_[pI] != -1)
            {
                data[pointMap_[pI]] = pointData[pI];
            }
        }
    }
}

void Foam::simpleVTKWriter::addPointData
(
    const word &name,
    const List<vector> &pointData
)
{
    if (canAddPointData_)
    {
        pointVectorNames_.append(name);
        pointVectorData_.append(List<vector>(points_.size()));

        List<vector> &data = pointVectorData_[pointVectorData_.size()-1];
        forAll(pointData,pI)
        {
            if (pointMap_[pI] != -1)
            {
                data[pointMap_[pI]] = pointData[pI];
            }
        }
    }
}

void Foam::simpleVTKWriter::addPointDataToPoints
(
    const word &name,
    const List<label> &pointData
)
{
    pointLabelNames_.append(name);
    pointLabelData_.append(List<label>(points_.size()));

    pointLabelData_[pointLabelData_.size()-1] = pointData;
}


void Foam::simpleVTKWriter::addPointDataToPoints
(
    const word &name,
    const List<scalar> &pointData
)
{
    pointScalarNames_.append(name);
    pointScalarData_.append(List<scalar>(points_.size()));

    pointScalarData_[pointScalarData_.size()-1] = pointData;
}

void Foam::simpleVTKWriter::addPointDataToPoints
(
    const word &name,
    const List<vector> &pointData
)
{
    pointVectorNames_.append(name);
    pointVectorData_.append(List<vector>(points_.size()));

    pointVectorData_[pointVectorData_.size()-1] = pointData;
}


void Foam::simpleVTKWriter::write
(
    const fileName &fileName,
    const bool parData
) const
{
    OFstream out(fileName);

    if (Pstream::parRun() && parData)
    {
        // Gather all points
        List<pointField> gatheredPoints(Pstream::nProcs());
        gatheredPoints[Pstream::myProcNo()] = points_;
        Pstream::gatherList(gatheredPoints);

        // Gather all faces
        List<List<face>> gatheredFaces(Pstream::nProcs());
        gatheredFaces[Pstream::myProcNo()] = faces_;
        Pstream::gatherList(gatheredFaces);

        List<point> globalPoints(0);
        List<face> globalFaces(0);

        if (Pstream::master())
        {
            // On master combine all points
            globalPoints =
                ListListOps::combine<pointField>
                (
                    gatheredPoints,
                    accessOp<pointField>()
                );

            label pointOffset = 0;
            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
                List<face>& localFaces = gatheredFaces[procI];
                forAll(localFaces, faceI)
                {
                    face& f = localFaces[faceI];
                    forAll(f, fp)
                    {
                        f[fp] += pointOffset;
                    }
                }
                pointOffset += gatheredPoints[procI].size();
            }

            globalFaces =
                ListListOps::combine<List<face>>
                (
                    gatheredFaces,
                    accessOp<List<face>>()
                );

            out << "# vtk DataFile Version 2.0\n";
            out << "Surface\n";
            out << "ASCII\n";
            out << "DATASET POLYDATA\n";
            out << "POINTS " << globalPoints.size() << " float\n";

            forAll(globalPoints,pI)
            {
                out << globalPoints[pI][0] << " " << globalPoints[pI][1] << " "
                    << globalPoints[pI][2] << "\n";
            }

            out << "POLYGONS ";

            label nFaceNodes = 0;
            forAll(globalFaces,fI)
            {
                nFaceNodes += globalFaces[fI].size();
            }

            out << globalFaces.size() << " " << globalFaces.size() + nFaceNodes
                << "\n";

            forAll(globalFaces,fI)
            {
                out << globalFaces[fI].size();
                forAll(globalFaces[fI],pI)
                {
                    out << " " << globalFaces[fI][pI];
                }
                out << "\n";
            }
        }

        if (Pstream::master())
        {
            if
            (
                !faceLabelData_.empty() || !faceScalarData_.empty()
                || !faceVectorData_.empty()
            )
            {
                out << "CELL_DATA " << globalFaces.size() << "\n";
                out << "FIELD attributes "
                    << faceLabelData_.size() + faceScalarData_.size()
                    + faceVectorData_.size() << "\n";
            }
        }

        forAll(faceLabelData_,setI)
        {
            List<Field<label>> gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = faceLabelData_[setI];
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<label> allValues
                (
                    ListListOps::combine<Field<label>>
                    (
                        gatheredValues,
                        accessOp<Field<label>>()
                     )
                 );
                out << faceLabelNames_[setI] << " 1 "
                    << allValues.size()
                    << " float\n"; //this could be made more flexible
                forAll(allValues,i)
                {
                    out << allValues[i] << "\n";
                }
            }
        }

        forAll(faceScalarData_,setI)
        {
            List<Field<scalar>> gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = faceScalarData_[setI];
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<scalar> allValues
                (
                    ListListOps::combine<Field<scalar>>
                    (
                        gatheredValues,
                        accessOp<Field<scalar>>()
                     )
                 );
                out << faceScalarNames_[setI] << " 1 "
                    << allValues.size()
                    << " float\n"; //this could be made more flexible
                forAll(allValues,i)
                {
                    out << allValues[i] << "\n";
                }
            }
        }

        forAll(faceVectorData_,setI)
        {
            List<Field<vector>> gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = faceVectorData_[setI];
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<vector> allValues
                (
                    ListListOps::combine<Field<vector>>
                    (
                        gatheredValues,
                        accessOp<Field<vector>>()
                     )
                 );
                out << faceVectorNames_[setI] << " 3 "
                    << allValues.size()
                    << " float\n"; //this could be made more flexible
                forAll(allValues,i)
                {
                    out << allValues[i][0] << " "
                        << allValues[i][1] << " "
                        << allValues[i][2] << "\n";
                }
            }
        }
        if (Pstream::master())
        {
            if
            (
                !pointLabelData_.empty() || !pointScalarData_.empty()
                || !pointVectorData_.empty()
            )
            {
                out << "POINT_DATA " << globalPoints.size() << "\n";
                out << "FIELD attributes "
                    << pointLabelData_.size() + pointScalarData_.size()
                    + pointVectorData_.size() << "\n";
            }
        }

        forAll(pointLabelData_,setI)
        {
            List<Field<label>> gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = pointLabelData_[setI];
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<label> allValues
                (
                    ListListOps::combine<Field<label>>
                    (
                        gatheredValues,
                        accessOp<Field<label>>()
                     )
                 );
                out << pointLabelNames_[setI] << " 1 "
                    << allValues.size()
                    << " float\n"; //this could be made more flexible
                forAll(allValues,i)
                {
                    out << allValues[i] << "\n";
                }
            }
        }

        forAll(pointScalarData_,setI)
        {
            List<Field<scalar>> gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = pointScalarData_[setI];
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<scalar> allValues
                (
                    ListListOps::combine<Field<scalar>>
                    (
                        gatheredValues,
                        accessOp<Field<scalar>>()
                     )
                 );
                out << pointScalarNames_[setI] << " 1 "
                    << allValues.size()
                    << " float\n"; //this could be made more flexible
                forAll(allValues,i)
                {
                    out << allValues[i] << "\n";
                }
            }
        }

        forAll(pointVectorData_,setI)
        {
            List<Field<vector>> gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = pointVectorData_[setI];
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<vector> allValues
                (
                    ListListOps::combine<Field<vector>>
                    (
                        gatheredValues,
                        accessOp<Field<vector>>()
                     )
                 );
                out << pointVectorNames_[setI] << " 3 "
                    << allValues.size()
                    << " float\n"; //this could be made more flexible
                forAll(allValues,i)
                {
                    out << allValues[i][0] << " "
                        << allValues[i][1] << " "
                        << allValues[i][2] << "\n";
                }
            }
        }
    }
    else
    {

        out << "# vtk DataFile Version 2.0\n";
        out << "Surface\n";
        out << "ASCII\n";
        out << "DATASET POLYDATA\n";
        out << "POINTS " << points_.size() << " float\n";

        forAll(points_,pI)
        {
            out << points_[pI][0] << " " << points_[pI][1] << " "
                << points_[pI][2] << "\n";
        }

        out << "POLYGONS ";

        label nFaceNodes = 0;
        forAll(faces_,fI)
        {
            nFaceNodes += faces_[fI].size();
        }

        out << faces_.size() << " " << faces_.size() + nFaceNodes << "\n";

        forAll(faces_,fI)
        {
            out << faces_[fI].size();
            forAll(faces_[fI],pI)
            {
                out << " " << faces_[fI][pI];
            }
            out << "\n";
        }

        if
        (
            !faceLabelData_.empty() || !faceScalarData_.empty()
            || !faceVectorData_.empty()
        )
        {
            out << "CELL_DATA " << faces_.size() << "\n";
            out << "FIELD attributes "
                << faceLabelData_.size() + faceScalarData_.size()
                + faceVectorData_.size() << "\n";

            forAll(faceLabelData_,setI)
            {
                out << faceLabelNames_[setI] << " 1 "
                    << faceLabelData_[setI].size()
                    << " float\n"; //this could be made more flexible
                forAll(faceLabelData_[setI],i)
                {
                    out << faceLabelData_[setI][i] << "\n";
                }
            }

            forAll(faceScalarData_,setI)
            {
                out << faceScalarNames_[setI] << " 1 "
                    << faceScalarData_[setI].size()
                    << " float\n"; //this could be made more flexible
                forAll(faceScalarData_[setI],i)
                {
                    out << faceScalarData_[setI][i] << "\n";
                }
            }

            forAll(faceVectorData_,setI)
            {
                //this could be made more flexible
                out << faceVectorNames_[setI] << " 3 "
                    << faceVectorData_[setI].size() << " float\n";
                forAll(faceVectorData_[setI],i)
                {
                    out << faceVectorData_[setI][i][0] << " "
                        << faceVectorData_[setI][i][1] << " "
                        << faceVectorData_[setI][i][2] << "\n";
                }
            }
        }

        if
        (
            !pointLabelData_.empty() || !pointScalarData_.empty()
            || !pointVectorData_.empty()
        )
        {

            out << "POINT_DATA " << points_.size() << "\n";
            out << "FIELD attributes "
                << pointLabelData_.size() + pointScalarData_.size()
                + pointVectorData_.size() << "\n";

            forAll(pointLabelData_,setI)
            {
                //this could be made more flexible
                out << pointLabelNames_[setI] << " 1 "
                    << pointLabelData_[setI].size() << " float\n";
                forAll(pointLabelData_[setI],i)
                {
                    out << pointLabelData_[setI][i] << "\n";
                }
            }
            forAll(pointScalarData_,setI)
            {
                //this could be made more flexible
                out << pointScalarNames_[setI] << " 1 "
                    << pointScalarData_[setI].size() << " float\n";
                forAll(pointScalarData_[setI],i)
                {
                    out << pointScalarData_[setI][i] << "\n";
                }
            }
            forAll(pointVectorData_,setI)
            {
                    //this could be made more flexible
                out <<  pointVectorNames_[setI] << " 3 "
                    << pointVectorData_[setI].size() << " float\n";
                forAll(pointVectorData_[setI],i)
                {
                    out << pointVectorData_[setI][i][0] << " "
                        << pointVectorData_[setI][i][1] << " "
                        << pointVectorData_[setI][i][2] << "\n";
                }
            }
        }
    }
}
