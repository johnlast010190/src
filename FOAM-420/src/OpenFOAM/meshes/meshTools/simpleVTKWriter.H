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

Class
    Foam::simpleVTKWriter

Description
    Dumps face lists out as VTK file. Intended for visual debugging of mesh
    generation. As default, eliminates nodes that aren't referenced by any
    faces. Accepts face and point scalar and vector data. Point data is assumed
    to be supplied on the same numbering basis as the original points.

SourceFiles
    simpleVTKWriter.C

\*---------------------------------------------------------------------------*/

#ifndef simpleVTKWriter_H
#define simpleVTKWriter_H

#include "containers/Lists/List/List.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/meshShapes/face/face.H"
#include "primitives/strings/fileName/fileName.H"

namespace Foam
{
    class boundBox;

    class simpleVTKWriter
    {
        private:

        List<face> faces_;
        pointField points_;
        labelList pointMap_;
        labelList faceMap_;

        const bool canAddPointData_;
        const bool canAddFaceData_;

        DynamicList<List<label>> faceLabelData_;
        DynamicList<List<scalar>> faceScalarData_;
        DynamicList<List<vector>> faceVectorData_;

        DynamicList<List<label>> pointLabelData_;
        DynamicList<List<scalar>> pointScalarData_;
        DynamicList<List<vector>> pointVectorData_;

        DynamicList<word> faceLabelNames_;
        DynamicList<word> faceScalarNames_;
        DynamicList<word> faceVectorNames_;

        DynamicList<word> pointLabelNames_;
        DynamicList<word> pointScalarNames_;
        DynamicList<word> pointVectorNames_;


        public:

        simpleVTKWriter();

        simpleVTKWriter
        (
            const List<face> &faces,
            const pointField &points,
            const bool removeUnusedPoints = true
        );

        simpleVTKWriter
        (
            const List<face> &faces,
            const labelList &faceList,
            const pointField &points
        );

        simpleVTKWriter(const labelList &pointSet, const List<vector> &points);

        simpleVTKWriter
        (
            const List<labelList> &pointSet,
            const List<vector> &points
         );

        simpleVTKWriter(const List<vector> &points);

        simpleVTKWriter(const List<edge> &edges, const List<vector> &points);

        simpleVTKWriter
        (
            const List<edge> &edges,
            const List<vector> &points,
            const DynamicList<label> &edgeList
         );

        simpleVTKWriter(const boundBox &bb);

        void addMorePoints(const List<vector> &points);

        void addMorePoints
        (
            const labelList &pointList,
            const List<vector> &points
         );


        void addFaceData(const word &name, const labelList &faceData);
        void addFaceData(const word &name, const List<scalar> &faceData);
        void addFaceData(const word &name, const List<vector> &faceData);

        void addPointData(const word &name, const labelList &pointData);
        void addPointData(const word &name, const List<scalar> &pointData);
        void addPointData(const word &name, const List<vector> &pointData);

        void addPointDataToPoints(const word &name, const labelList &pointData);
        void addPointDataToPoints(const word &name, const List<scalar> &pointData);
        void addPointDataToPoints(const word &name, const List<vector> &pointData);

        void write
        (
            const fileName &filename,
            const bool parData = true
        ) const;
    };

}

#endif
