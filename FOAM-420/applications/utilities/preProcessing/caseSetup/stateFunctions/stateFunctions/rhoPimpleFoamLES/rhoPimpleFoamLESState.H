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
    Foam::rhoPimpleFoamLESState

Description
    pimpleFoam state identifiers and specialised assembly functions

SourceFiles
    rhoPimpleFoamLESState.C

\*---------------------------------------------------------------------------*/

#ifndef rhoPimpleFoamLESState_H
#define rhoPimpleFoamLESState_H

#include "stateFunctions/rhoPimpleFoamRAS/rhoPimpleFoamRASState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rhoPimpleFoamLESState Declaration
\*---------------------------------------------------------------------------*/

class rhoPimpleFoamLESState
:
    public rhoPimpleFoamRASState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        rhoPimpleFoamLESState(const rhoPimpleFoamLESState&);

        //- Disallow default bitwise assignment
        void operator=(const rhoPimpleFoamLESState&);


public:

    //- Runtime type information
        TypeName("rhoPimpleFoamLES");

    // Constructors

        //- Construct from components
        rhoPimpleFoamLESState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~rhoPimpleFoamLESState()
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
