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
    (c) 2016 Esi Ltd.

Class
    Foam::pimpleFoamLESState

Description
    pimpleFoam state identifiers and specialised assembly functions

SourceFiles
    pimpleFoamLESState.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleFoamLESState_H
#define pimpleFoamLESState_H

#include "stateFunctions/pimpleFoamRAS/pimpleFoamRASState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class pimpleFoamLESState Declaration
\*---------------------------------------------------------------------------*/

class pimpleFoamLESState
:
    public pimpleFoamRASState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        pimpleFoamLESState(const pimpleFoamLESState&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleFoamLESState&);


public:

    //- Runtime type information
        TypeName("pimpleFoamLES");

    // Constructors

        //- Construct from components
        pimpleFoamLESState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~pimpleFoamLESState()
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
