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
    (c) 2019-2022 Esi Ltd.

Class
    Foam::trapDoor

Description
    Active patch based on pressure drop
    Deactivate patch based on normal velocity of fluid
    If the alpha field does not exist, take the velocity of the whole patch
    into accout which will work for single phase flow

SourceFiles
    trapDoor.C

\*---------------------------------------------------------------------------*/

#ifndef trapDoor_H
#define trapDoor_H

#include "patchSwitchFvMesh/patchSwitches/patchSwitch/patchSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class trapDoor Declaration
\*---------------------------------------------------------------------------*/

class trapDoor
:
    public patchSwitch
{
    // Private Member Functions

        //- Direction invertion
        Switch invertDirection_;

        //- Velocity Threashold
        scalar fVel_;

        //- Pressure drop Threashold
        scalar dp_;

        //- Name of fraction field
        word UName_;

        //- Name of fraction field
        word pName_;

        //- Name of fraction field
        word fractionName_;

        virtual bool enableCondition();

        virtual bool disableCondition();

        //- Disallow default bitwise copy construct
        trapDoor(const trapDoor&);

        //- Disallow default bitwise assignment
        void operator=(const trapDoor&);


public:

    //- Runtime type information
    TypeName("trapDoor");


    // Constructors

        //- Construct from objectRegistry, and read/write options
        trapDoor
        (
            const fvMesh&,
            const dictionary&
        );

    // Selectors


    //- Destructor
    virtual ~trapDoor(){};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
