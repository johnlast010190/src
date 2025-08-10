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
    (c) 2018 Esi Ltd.

Class
    Foam::constraintBoundaryPoints

Description
   Constrain the patch points to move on the edge direction

SourceFiles
    constraintBoundaryPoints.C

\*---------------------------------------------------------------------------*/

#ifndef constraintBoundaryPoints_H
#define constraintBoundaryPoints_H

#include "twoDPointCorrector/twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class constraintBoundaryPoints Declaration
\*---------------------------------------------------------------------------*/

class constraintBoundaryPoints
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& pMesh_;

        //- Reference to basePoints
        const pointField& basePoints_;

        //- Reference to basePoints
        const vectorField& baseSf_;

        //- patch points
        boolList boundaryPoints_;

        //- gib points
        boolList gibPoints_;

        //- gib && boundary points
        boolList gibBoundaryPoints_;


    // Private functions

        //-
        void createDataStructure();

        //-
        void markBoundaryPoints();

        //-
        void markGIBPoints(const labelList&);

public:

    // Constructors

        //- Construct from components
        constraintBoundaryPoints
        (
            const polyMesh&,
            const pointField&,
            const vectorField&
        );

    //- Destructor
    ~constraintBoundaryPoints(){};

    // Member Functions

        //-
        void correctPoints
        (
            pointField&,
            const labelList&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
