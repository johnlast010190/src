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
    Allow calculation and storage of global edge cells

SourceFiles
    extendedEdgeCells.C

\*---------------------------------------------------------------------------*/

#ifndef extendedEdgeCells_H
#define extendedEdgeCells_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshes/primitiveShapes/point/pointField.H"
#include "containers/HashTables/HashSet/HashSet.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class extendedEdgeCells Declaration
\*---------------------------------------------------------------------------*/


class extendedEdgeCells
{
    //- List of edgeCells in global addressing
    labelList edgeCells_;

    //- List of edge cell centres
    pointField edgeCellCentres_;

    //- List of edgeFaces face owner and neighbour in global addressing
    List<labelPair> edgeFacesOwnNei_;

    //- Global edge ID
    label globalEdgeID_;

    public:

    class sumEqOp
    {
        public:
        void operator()
        ( extendedEdgeCells& x, const extendedEdgeCells& y ) const
        {
            if (x.edgeCells().size() == 0)
            {
                x.setCells(y.edgeCells());
                x.setCentres(y.edgeCentres());
                x.setOwnNei(y.edgeFacesOwnNei());
            }
            else
            {
                labelList& mEdgeCells = x.edgeCells();
                pointField& mEdgeCentres = x.edgeCentres();
                List<labelPair>& mEdgeFacesOwnNei = x.edgeFacesOwnNei();
                const label mGEID = x.globalEdgeID();

                const labelList& sEdgeCells = y.edgeCells();
                const pointField& sEdgeCentres = y.edgeCentres();
                const List<labelPair>& sEdgeFacesOwnNei = y.edgeFacesOwnNei();
                const label sGEID = y.globalEdgeID();

                if (mGEID != sGEID)
                {
                    label sz = mEdgeCells.size();
                    label newSz = sz + sEdgeCells.size();
                    mEdgeCells.setSize(newSz);
                    mEdgeCentres.setSize(newSz);

                    label j = 0;
                    for (label i = sz; i < newSz; i++)
                    {
                        mEdgeCells[i] = sEdgeCells[j];
                        mEdgeCentres[i] = sEdgeCentres[j];
                        j++;
                    }

                    sz = mEdgeFacesOwnNei.size();
                    newSz = sz + sEdgeFacesOwnNei.size();
                    mEdgeFacesOwnNei.setSize(newSz);

                    j = 0;
                    for (label i = sz; i < newSz; i++)
                    {
                        mEdgeFacesOwnNei[i] = sEdgeFacesOwnNei[j];
                        j++;
                    }
                }
            }
        }
    };

    // Construct from components
    extendedEdgeCells
    (
        labelList ec=labelList(0),
        pointField ecc=pointField(0),
        List<labelPair> efon=List<labelPair>(0),
        label geID=-1
    )
    :
        edgeCells_(ec),
        edgeCellCentres_(ecc),
        edgeFacesOwnNei_(efon),
        globalEdgeID_(geID)
    {};

    // Access functions
    labelList& edgeCells() {return edgeCells_;};
    pointField& edgeCentres() {return edgeCellCentres_;};
    List<labelPair>& edgeFacesOwnNei() {return edgeFacesOwnNei_;};

    labelList edgeCells() const {return edgeCells_;};
    pointField edgeCentres() const {return edgeCellCentres_;};
    List<labelPair> edgeFacesOwnNei() const {return edgeFacesOwnNei_;};
    label globalEdgeID() const {return globalEdgeID_;};

    // Manipulation functions
    void setCells(labelList ec){edgeCells_=ec;}
    void setCentres(pointField ecc){edgeCellCentres_=ecc;}
    void setOwnNei(List<labelPair> efon){edgeFacesOwnNei_=efon;}

    //- Construct from Istream.
    extendedEdgeCells(Istream&);

    friend Istream& operator>>
    (Istream& is, extendedEdgeCells&);
    friend Ostream& operator<<
    (Ostream& os, const extendedEdgeCells&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

