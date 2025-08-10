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
    (c) 2014 Esi Ltd.

Class
    Foam::activeSet

Description
    Class that keeps information for the active sets of points, cells, faces and
    edges that are currently being optimized.

SourceFiles
    activeSet.C

\*---------------------------------------------------------------------------*/

#ifndef activeSet_H
#define activeSet_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class activeSet Declaration
\*---------------------------------------------------------------------------*/

class activeSet
{
    // Private Data

        labelList pointSet_;
        labelList faceSet_;
        labelList cellSet_;
        labelList faceCellSet_;
        labelList edgeSet_;
        bool updated_;

public    :
        //Construct given the fields
        activeSet
        (
            const labelList& pointSet,
            const labelList& faceSet,
            const labelList& cellSet,
            const labelList& faceCellSet,
            const labelList& edgeSet
        );

        //Construct only from points
        activeSet
        (
            const labelList& pointSet
        );

        //construct null
        activeSet();

    virtual ~activeSet();

    //Member Function

    inline const labelList& activePoints() const
    {
        return pointSet_;
    };

    inline const labelList& activeFaces() const
    {
        return faceSet_;
    };

    inline const labelList& activeFaceCells() const
    {
        return faceCellSet_;
    };

    inline const labelList& activeCells() const
    {
        return cellSet_;
    };

    inline const labelList& activeEdges() const
    {
        return edgeSet_;
    };

    void update
    (
        const labelList& pointSet,
        const labelList& faceSet,
        const labelList& cellSet,
        const labelList& faceCellSet,
        const labelList& edgeSet
    );

    void updatePointSet(const labelList& pointSet);

    bool updated() const;
};

} /* namespace Foam */

#endif /* activeSet_H */
