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
    (c) 2019 Esi Ltd.

Class
    Foam::tableSwitch

Description
    Activate/Deactivate based on userDefined table

SourceFiles
    tableSwitch.C

\*---------------------------------------------------------------------------*/

#ifndef tableSwitch_H
#define tableSwitch_H

#include "patchSwitchFvMesh/patchSwitches/patchSwitch/patchSwitch.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tableSwitch Declaration
\*---------------------------------------------------------------------------*/

class tableSwitch
:
    public patchSwitch
{

        autoPtr<Function1<bool>> boolTable_;

    // Private Member Functions

        virtual bool enableCondition();

        virtual bool disableCondition();

        //- Disallow default bitwise copy construct
        tableSwitch(const tableSwitch&);

        //- Disallow default bitwise assignment
        void operator=(const tableSwitch&);


public:

    //- Runtime type information
    TypeName("tableSwitch");


    // Constructors

        //- Construct from objectRegistry, and read/write options
        tableSwitch
        (
            const fvMesh&,
            const dictionary&
        );

    // Selectors


    //- Destructor
    virtual ~tableSwitch(){};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
