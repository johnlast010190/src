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

Description

\*---------------------------------------------------------------------------*/

#ifndef hexBlock_H
#define hexBlock_H

#include "primitives/ints/lists/labelList.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/meshShapes/face/faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hexBlock Declaration
\*---------------------------------------------------------------------------*/

class hexBlock
{
    // Private data

        //- Handedness of the block
        enum handed
        {
            noPoints,
            right,
            left
        };

        //- Number of point in each direction
        label xDim_;
        label yDim_;
        label zDim_;

        //- Handedness of the block
        handed blockHandedness_;

        //- List of points
        pointField points_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        hexBlock(const hexBlock&);

        //- Disallow default bitwise assignment
        void operator=(const hexBlock&);

        //- Vertex addressing inside the block
        inline label vtxLabel(label i, label j, label k) const;

        //- Calculate handedness of block
        void setHandedness();

public:

    // Constructors

        //- Construct from components
        hexBlock
        (
            const label nx,
            const label ny,
            const label nz
        );

    // Member Functions

        //- Number of points
        label xDim() const
        {
            return xDim_;
        }

        label yDim() const
        {
            return yDim_;
        }

        label zDim() const
        {
            return zDim_;
        }

        label nBlockPoints() const
        {
            return (xDim_ + 1)*(yDim_ + 1)*(zDim_ + 1);
        }

        label nBlockCells() const
        {
            return xDim_*yDim_*zDim_;
        }

        //- Return block points
        const pointField& points() const
        {
            if (blockHandedness_ == noPoints)
            {
                FatalErrorInFunction
                    << "points not read in yet"
                    << abort(FatalError);
            }

            return points_;
        }

        //- Return block cells
        labelListList blockCells() const;

        //- Return block patch faces given direction and range limits
        // From the cfx manual: direction
        // 0 = solid (3-D patch),
        // 1 = high i, 2 = high j, 3 = high k
        // 4 = low i, 5 = low j, 6 = low k
        faceList patchFaces(label direc, const labelList& range) const;


        //- Read block points either with or without blanking after every block.
        //  If twoDThickness > 0 reads (half) the points and extrudes the
        //  points in z direction.
        void readPoints
        (
            const bool readBlank,
            const scalar twoDThicknes,
            Istream&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
