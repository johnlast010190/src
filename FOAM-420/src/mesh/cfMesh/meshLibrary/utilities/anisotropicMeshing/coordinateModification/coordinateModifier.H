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
    coordinateModifier

Description
    Base class for modification of point coordinates
    in surface and volume meshes.

SourceFiles
    coordinateModifier.C

\*---------------------------------------------------------------------------*/

#ifndef coordinateModifier_H
#define coordinateModifier_H

#include "primitives/strings/word/word.H"
#include "meshes/primitiveShapes/point/point.H"
#include "utilities/anisotropicMeshing/coordinateModification/coordinateModificationList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class boundBox;

/*---------------------------------------------------------------------------*\
                   Class coordinateModifier Declaration
\*---------------------------------------------------------------------------*/

class coordinateModifier
{
    // Private data

        //- Reference to dictionary
        const dictionary& modificationDict_;

        //- list of coordinate modification objects
        coordinateModificationList modifiers_;

        //- list of backward coordinate modifiers
        coordinateModificationList backwardModifiers_;

    // Private member functions
        //- check if the comibination of modifiers has a valid inverse
        void checkForValidInverse() const;

public:

    // Constructors

        //- Construct from dictionary
        coordinateModifier(const dictionary&);

    // Destructor

        ~coordinateModifier();


    // Member Functions

        //- calculate the modified coordinate of the point
        point modifiedPoint(const point&) const;

        //- calculate the displacement vector for the backward modification
        point backwardModifiedPoint(const point&) const;

        //- Print modification objects
        void printObjects() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
