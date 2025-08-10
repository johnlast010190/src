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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::polyModifyPoint

Description
    Class describing modification of a point.

\*---------------------------------------------------------------------------*/

#ifndef polyModifyPoint_H
#define polyModifyPoint_H

#include "primitives/ints/label/label.H"
#include "meshes/primitiveShapes/point/point.H"
#include "polyTopoChange/polyTopoChange/topoAction/topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class polyModifyPoint Declaration
\*---------------------------------------------------------------------------*/

class polyModifyPoint
:
    public topoAction
{
    // Private data

        //- Point ID
        label pointID_;

        //- New point location
        point location_;

        //- Remove from current zone
        bool removeFromZone_;

        //- New zone ID
        label zoneID_;

        //- Does the point support a cell
        bool inCell_;


public:

    // Static data members

        //- Runtime type information
        TypeName("modifyPoint");


    // Constructors

        //- Construct null.  Used only for list construction
        polyModifyPoint()
        :
            pointID_(-1),
            location_(Zero),
            removeFromZone_(false),
            zoneID_(-1),
            inCell_(false)
        {}

        //- Construct from components
        polyModifyPoint
        (
            const label pointID,
            const point& newP,
            const bool removeFromZone,
            const label newZoneID,
            const bool inCell
        )
        :
            pointID_(pointID),
            location_(newP),
            removeFromZone_(removeFromZone),
            zoneID_(newZoneID),
            inCell_(inCell)
        {}

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyModifyPoint(*this));
        }


    // Default Destructor

    // Member Functions

        //- Point ID
        label pointID() const
        {
            return pointID_;
        }

        //- New point location
        const point& newPoint() const
        {
            return location_;
        }

        //- Does the point belong to a zone?
        bool isInZone() const
        {
            return zoneID_ >= 0;
        }

        //- Should the point be removed from current zone
        bool removeFromZone() const
        {
            return removeFromZone_;
        }

        //- Point zone ID
        label zoneID() const
        {
            return zoneID_;
        }

        //- Does the point support a cell
        bool inCell() const
        {
            return inCell_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
