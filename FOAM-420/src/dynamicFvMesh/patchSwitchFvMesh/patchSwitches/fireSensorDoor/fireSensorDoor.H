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
    Foam::fireSensorDoor

Description
    Avtivate GIB patch if the the scalarField exceeds a value at a face of
    patched.
    It can be extended the sensor to look values of a cellZone etc etc

    If it GIB is activated is not getting disabled

SourceFiles
    fireSensorDoor.C

\*---------------------------------------------------------------------------*/

#ifndef fireSensorDoor_H
#define fireSensorDoor_H

#include "patchSwitchFvMesh/patchSwitches/patchSwitch/patchSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fireSensorDoor Declaration
\*---------------------------------------------------------------------------*/

class fireSensorDoor
:
    public patchSwitch
{
    // Private Member Functions

        //- Value of threshold
        const scalar threshold_;

        //- Patches to investigate
        wordList patchNames_;

        //- Points to investigate
        List<point> points_;

        //- Name of field
        word psName_;

        virtual bool enableCondition();

        virtual bool disableCondition();

        //- Disallow default bitwise copy construct
        fireSensorDoor(const fireSensorDoor&);

        //- Disallow default bitwise assignment
        void operator=(const fireSensorDoor&);


public:

    //- Runtime type information
    TypeName("fireSensorDoor");


    // Constructors

        //- Construct from objectRegistry, and read/write options
        fireSensorDoor
        (
            const fvMesh&,
            const dictionary&
        );

    // Selectors


    //- Destructor
    virtual ~fireSensorDoor(){};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
