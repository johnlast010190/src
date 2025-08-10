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
    (c) 2016-2017 Esi Ltd.

Class
    Foam::autoHexDualiser

Description
    Creates dual of polyMesh. No boundary handling

SourceFiles
    autoHexDualiser.C

\*---------------------------------------------------------------------------*/

#ifndef autoHexDualiser_H
#define autoHexDualiser_H

#include "meshRefinement/meshRefinement.H"
#include "autoHexDualiser/extendedEdgeCells.H"
#include "autoHexDualiser/extendedPointEdgeFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoHexDualiser Declaration
\*---------------------------------------------------------------------------*/

class autoHexDualiser
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reconstruct edge faces
        void calcFaceLoop
        (
            const extendedEdgeCells& eEC,
            face& f,
            pointField& pts
        );

        //- Calculate boundary zones after dualisation
        void calcBoundaryZones
        (
            labelList& pointCellZones
        );

        //- Create dual mesh
        void calcDual();

        //Update gap cell field if present
        void updateGapCells();

        //- post splitting boundary faces to avoid degenerate cells
        void splitBoundaryFaces();

        class unionEqOp
        {
        public:
            void operator()
            (
                List<labelPair>& x,
                const List<labelPair>& y
            ) const
            {
                if (y.size())
                {
                    if (x.empty())
                    {
                        x = y;
                    }
                    else
                    {
                        x.append(y);
                    }
                }
            }
        };


public:

    //- Runtime type information
    ClassName("autoHexDualiser");


    // Constructors
    autoHexDualiser(meshRefinement&);

    //- Destructor
    ~autoHexDualiser();


    // Member Functions

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
