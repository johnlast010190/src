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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2012 OpenFOAM Foundation

Class
    Foam::tetPoints

Description
    Tet storage. Null constructable (unfortunately tetrahedron<point, point>
    is not)

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef tetPoints_H
#define tetPoints_H

#include "meshes/primitiveShapes/tetrahedron/tetrahedron.H"
#include "containers/Lists/FixedList/FixedList.H"
#include "meshes/treeBoundBox/treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tetPoints Declaration
\*---------------------------------------------------------------------------*/

class tetPoints
:
    public FixedList<point, 4>
{
public:

    // Constructors

        //- Construct null
        inline tetPoints()
        {}

        //- Construct from four points
        inline tetPoints
        (
            const point& a,
            const point& b,
            const point& c,
            const point& d
        )
        {
            operator[](0) = a;
            operator[](1) = b;
            operator[](2) = c;
            operator[](3) = d;
        }

    // Member Functions

        //- Return the tetrahedron
        inline tetPointRef tet() const
        {
            return tetPointRef
            (
                operator[](0),
                operator[](1),
                operator[](2),
                operator[](3)
            );
        }

        //- Calculate the bounding box
        inline treeBoundBox bounds() const
        {
            treeBoundBox bb(operator[](0));
            for (label i = 1; i < size(); ++i)
            {
                bb.add(operator[](i));
            }
            return bb;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
