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

Description
    Class for checking if points are crossing processor boundaries and they are
    hiting the boundaries. If a pointHit for a point is found then snap on the
    boundary

SourceFiles
    parallelIntersectionData.C

\*---------------------------------------------------------------------------*/

#ifndef parallelIntersectionData_H
#define parallelIntersectionData_H

#include "db/typeInfo/typeInfo.H"
#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "dynamicGIBFvMesh/movingGIBTools/intersectionData/intersectionData.H"
#include "containers/LinkedLists/user/LIFOStack.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class parallelIntersectionData Declaration
\*---------------------------------------------------------------------------*/

class parallelIntersectionData
{
    // Private data

    // reference to mesh
    const dynamicGIBFvMesh& mesh_;

    // base mesh points
    const pointField& basePoints_;

    // base mesh face centers
    const vectorField& baseCf_;

    // problematic (crossing boundaries) point data in each processor
    const List<intersectionData>& iDl_;

    // reconstructed list which gets updated until all points investigated
    List<List<intersectionData>> giDl_;


    // Private Member Functions

    // intersection checking algorithm
    // face --> cell --> face .. walk until it hit boundary
    void faceCellFaceItersections
    (
        intersectionData&,
        LIFOStack<label>&,
        boolList&,
        boolList&
    );


    // prepare the intersection checking algorithm
    void itersectionChecking();

    // send problematic points
    void sendPointsToProcessors
    (
        DynamicList<intersectionData>&,
        DynamicList<label>&
    );

    // update points
    void receiveAndUpdateData
    (
        const List<intersectionData>&,
        const List<label>&
    );

    // check if there are any problematic points in every proc
    bool doLoop();


public:


    // Constructors

        //- Construct from components
        parallelIntersectionData
        (
            const dynamicGIBFvMesh&,
            const pointField&,
            const vectorField&,
            const List<intersectionData>&
        );


    //- Destructor
    virtual ~parallelIntersectionData(){};


    // Member Functions


        //- update the point displacement
        void modifyPointsPassingBoundary
        (
            pointField&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
