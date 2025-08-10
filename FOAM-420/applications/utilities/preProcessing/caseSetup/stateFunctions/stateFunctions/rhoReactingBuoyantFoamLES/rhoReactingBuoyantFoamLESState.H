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
    (c) 2016-2019 Esi Ltd.

Class
    Foam::rhoReactingBuoyantFoamLESState

Description
    pimpleFoam state identifiers and specialised assembly functions

SourceFiles
    rhoReactingBuoyantFoamLESState.C

\*---------------------------------------------------------------------------*/

#ifndef rhoReactingBuoyantFoamLESState_H
#define rhoReactingBuoyantFoamLESState_H

#include "stateFunctions/reactingFoamRAS/reactingFoamRASState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rhoReactingBuoyantFoamLESState Declaration
\*---------------------------------------------------------------------------*/

class rhoReactingBuoyantFoamLESState
:
    public reactingFoamRASState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        rhoReactingBuoyantFoamLESState(const rhoReactingBuoyantFoamLESState&);

        //- Disallow default bitwise assignment
        void operator=(const rhoReactingBuoyantFoamLESState&);


public:

    //- Runtime type information
        TypeName("rhoReactingBuoyantFoamLES");

    // Constructors

        //- Construct from components
        rhoReactingBuoyantFoamLESState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~rhoReactingBuoyantFoamLESState()
    {}


    // Member Functions

    // transient time type indicator
    inline virtual turbulenceType turbulence() const
    {
        return tuLES;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
