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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::polyRemoveFace

Description
    Class containing data for face removal.

\*---------------------------------------------------------------------------*/

#ifndef polyRemoveFace_H
#define polyRemoveFace_H

#include "primitives/ints/label/label.H"
#include "polyTopoChange/polyTopoChange/topoAction/topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class polyRemoveFace Declaration
\*---------------------------------------------------------------------------*/

class polyRemoveFace
:
    public topoAction
{
    // Private data

        //- Face ID
        label faceID_;

        //- Merge faceID or -1
        label mergeFaceID_;

public:

    // Static data members

        //- Runtime type information
        TypeName("removeFace");


    // Constructors

        //- Construct null.  Used for constructing lists
        polyRemoveFace()
        :
            faceID_(-1),
            mergeFaceID_(-1)
        {}

        //- Construct from components
        polyRemoveFace(const label faceID, const label mergeFaceID = -1)
        :
            faceID_(faceID),
            mergeFaceID_(mergeFaceID)
        {}

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyRemoveFace(*this));
        }


    // Default Destructor


    // Member Functions

        //- Return face ID
        label faceID() const
        {
            return faceID_;
        }

        //- Return merge face ID
        label mergeFaceID() const
        {
            return mergeFaceID_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
