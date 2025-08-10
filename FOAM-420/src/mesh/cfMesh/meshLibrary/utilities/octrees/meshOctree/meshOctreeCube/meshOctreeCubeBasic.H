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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    meshOctreeCubeBasic

Description
    A cube stores information needed for mesh generation

SourceFiles
    meshOctreeCubeBasic.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeCubeBasic_H
#define meshOctreeCubeBasic_H

#include "utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeCoordinates.H"
#include "primitives/contiguous/contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class meshOctreeCubeBasic Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeCubeBasic
: public meshOctreeCubeCoordinates
{
    // Private data
        //- cube type
        direction cubeType_;

        //- processor to which the cube belongs to
        short procNo_;

public:

    // Constructors
        //- Null constructor
        inline meshOctreeCubeBasic();

        //- Construct from coordinates
        explicit inline meshOctreeCubeBasic
        (
            const meshOctreeCubeCoordinates&
        );

        //- Construct from coordinates and cube type
        inline meshOctreeCubeBasic
        (
            const meshOctreeCubeCoordinates& cc,
            const direction cubeType,
            const short procNo = ALLPROCS
        );

    // Destructor

        ~meshOctreeCubeBasic()
        {}

    // Enumerators

        enum typesOfCubes
        {
            UNKNOWN = 1,
            OUTSIDE = 2,
            DATA = 4,
            INSIDE = 8,
            ALLPROCS = -2,
            OTHERPROC = -3
        };

    // Member functions
        //- return type
        inline direction cubeType() const;

        //- set cube type
        inline void setCubeType(const direction);

        //- return processor number
        inline short procNo() const;

        //- set processor number
        inline void setProcNo(const short);

        //- return coordinates in the octree
        inline const meshOctreeCubeCoordinates& coordinates() const;

        // Friend operators

        friend inline Ostream& operator<<
        (
            Ostream&,
            const meshOctreeCubeBasic&
        );

        friend inline Istream& operator>>
        (
            Istream&,
            meshOctreeCubeBasic&
        );
};

//- Specify data associated with meshOctreeCubeBasic type is contiguous
template<>
constexpr bool contiguous<meshOctreeCubeBasic>() {
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeBasicI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
