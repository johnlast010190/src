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
    (c) 2016-2023 Esi Ltd.

Class
    Foam::VOFState

Description
    Abstract, VOF state for USF - will eventually replace multiphaseVOF

SourceFiles
    VOFState.C

\*---------------------------------------------------------------------------*/

#ifndef VOFState_H
#define VOFState_H

#include "stateComponents/multiphaseMultiComponent/multiphaseMultiComponentState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class VOFState Declaration
\*---------------------------------------------------------------------------*/

class VOFState
:
    public multiphaseMultiComponentState
{

protected:

    // Protected Data

    // Protected member functions


private:

    // Private Data


    // Private Functions

        //- Disallow default bitwise copy construct
        VOFState(const VOFState&);

        //- Disallow default bitwise assignment
        void operator=(const VOFState&);



public:

    //- Runtime type information
        TypeName("VOF");

    // Constructors

        //- Construct from components
        VOFState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~VOFState()
    {}


    // Member Functions

        // Access

        // Check

        // Edit

            //- assemble field definitions
            virtual void initialise();

            //- assemble field definitions
            virtual void finalise();

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
