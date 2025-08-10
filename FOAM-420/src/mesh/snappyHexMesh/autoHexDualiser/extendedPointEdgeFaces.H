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
    (c) 2017 Esi Ltd.

Class
    Foam::extendedEdgeCells

Description
    Allow storage of cell faces in global storage

SourceFiles
    extendedEdgeCells.C

\*---------------------------------------------------------------------------*/

#ifndef extendedPointEdgeFaces_H
#define extendedPointEdgeFaces_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshes/primitiveShapes/point/pointList.H"
#include "meshes/meshShapes/face/faceList.H"
#include "containers/HashTables/HashSet/HashSet.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class extendedPointEdgeFaces Declaration
\*---------------------------------------------------------------------------*/

class extendedPointEdgeFaces
{
    //- List of edge faces in global addressing
    faceList edgeFaces_;

    //- List of points
    pointList pts_;

    //- List of points
    labelList addressing_;

    //- List of edge patch names
    wordList edgePatchNames_;

    //- List of edge neighbours
    labelList neighbours_;

    public:

    class sumEqOp
    {
        public:
        void operator()
        ( extendedPointEdgeFaces& x, const extendedPointEdgeFaces& y ) const
        {
            if (x.edgeFaces().size() == 0)
            {
                x.setFaces(y.edgeFaces());
                x.setCentres(y.points());
                x.setAddressing(y.addressing());
                x.setEdgeNames(y.edgePatchNames());
                x.setNeighbours(y.neighbours());
            }
            else
            {
                faceList& mEdgeFaces = x.edgeFaces();
                pointList& mPts = x.points();
                labelList& mAddressing = x.addressing();
                wordList& mPatchNames = x.edgePatchNames();
                labelList& mNeighbours = x.neighbours();

                const faceList& sEdgeFaces = y.edgeFaces();
                const pointList& sPts = y.points();
                const labelList& sAddressing = y.addressing();
                const wordList& sPatchNames = y.edgePatchNames();
                const labelList& sNeighbours = y.neighbours();

                label sz = mEdgeFaces.size();
                label newSz = sz + sEdgeFaces.size();
                mEdgeFaces.setSize(newSz);
                mPatchNames.setSize(newSz);
                mNeighbours.setSize(newSz);

                label j = 0;
                for (label i = sz; i < newSz; i++)
                {
                    mEdgeFaces[i] = sEdgeFaces[j];
                    mPatchNames[i] = sPatchNames[j];
                    mNeighbours[i] = sNeighbours[j];
                    j++;
                }

                sz = mPts.size();
                newSz = sz + sPts.size();
                mPts.setSize(newSz);
                mAddressing.setSize(newSz);
                j = 0;
                for (label i = sz; i < newSz; i++)
                {
                    mPts[i] = sPts[j];
                    mAddressing[i] = sAddressing[j];
                    j++;
                }

            }
        }
    };

    // Construct from components
    extendedPointEdgeFaces
    (
        faceList ef=faceList(0),
        pointList pts=pointList(0),
        labelList pa=labelList(0),
        wordList ep=wordList(0),
        labelList nei=labelList(0)
    )
    :
        edgeFaces_(ef),
        pts_(pts),
        addressing_(pa),
        edgePatchNames_(ep),
        neighbours_(nei)
    {};

    // Access functions
    faceList& edgeFaces() {return edgeFaces_;};
    pointList& points() {return pts_;};
    labelList& addressing() {return addressing_;};
    wordList& edgePatchNames() {return edgePatchNames_;};
    labelList& neighbours() {return neighbours_;};

    faceList edgeFaces() const {return edgeFaces_;};
    pointList points() const {return pts_;};
    labelList addressing() const {return addressing_;};
    wordList edgePatchNames() const {return edgePatchNames_;};
    labelList neighbours() const {return neighbours_;};

    // Manipulation functions
    void setFaces(faceList ef){edgeFaces_=ef;}
    void setCentres(pointList pts){pts_=pts;}
    void setAddressing(labelList addr){addressing_=addr;}
    void setEdgeNames(wordList en){edgePatchNames_=en;}
    void setNeighbours(labelList nei){neighbours_=nei;}

    //- Construct from Istream.
    extendedPointEdgeFaces(Istream&);


    friend Istream& operator>>
    (Istream& is, extendedPointEdgeFaces&);
    friend Ostream& operator<<
    (Ostream& os, const extendedPointEdgeFaces&);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
