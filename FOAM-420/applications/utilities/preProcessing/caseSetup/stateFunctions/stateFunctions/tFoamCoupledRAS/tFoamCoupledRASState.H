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
    (c) 2020 Esi Ltd.

Class
    Foam::tFoamCoupledRASState

Description
    foamCoupled state identifiers and specialised assembly functions

SourceFiles
    tFoamCoupledRASState.C

\*---------------------------------------------------------------------------*/

#ifndef tFoamCoupledRASState_H
#define tFoamCoupledRASState_H

#include "stateFunctions/foamCoupled/foamCoupledState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tFoamCoupledRASState Declaration
\*---------------------------------------------------------------------------*/

class tFoamCoupledRASState
:
    public foamCoupledState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        tFoamCoupledRASState(const tFoamCoupledRASState&);

        //- Disallow default bitwise assignment
        void operator=(const tFoamCoupledRASState&);



public:

    //- Runtime type information
        TypeName("tFoamCoupledRAS");

    // Constructors

        //- Construct from components
        tFoamCoupledRASState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~tFoamCoupledRASState()
    {}


    // Member Functions

        // Access

            // type indicators

            // time (transient/steady)
            inline virtual timeType time() const
            {
                return ttTrans;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
