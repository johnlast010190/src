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
    Foam::rhoPimpleFoamRASState

Description
    pimpleFoam state identifiers and specialised assembly functions

SourceFiles
    rhoPimpleFoamRASState.C

\*---------------------------------------------------------------------------*/

#ifndef rhoPimpleFoamRASState_H
#define rhoPimpleFoamRASState_H

#include "stateFunctions/rhoSimpleFoam/rhoSimpleFoamState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rhoPimpleFoamRASState Declaration
\*---------------------------------------------------------------------------*/

class rhoPimpleFoamRASState
:
    public rhoSimpleFoamState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        rhoPimpleFoamRASState(const rhoPimpleFoamRASState&);

        //- Disallow default bitwise assignment
        void operator=(const rhoPimpleFoamRASState&);


public:

    //- Runtime type information
        TypeName("rhoPimpleFoamRAS");

    // Constructors

        //- Construct from components
        rhoPimpleFoamRASState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~rhoPimpleFoamRASState()
    {}


    // Member Functions

    // transient time type indicator
    inline virtual timeType time() const
    {
        return ttTrans;
    }

    //- final manipulations
    virtual void finalise();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
