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
    Foam::pimpleFoamRASState

Description
    pimpleFoam state identifiers and specialised assembly functions

SourceFiles
    pimpleFoamRASState.C

\*---------------------------------------------------------------------------*/

#ifndef pimpleFoamRASState_H
#define pimpleFoamRASState_H

#include "stateFunctions/simpleFoam/simpleFoamState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class pimpleFoamRASState Declaration
\*---------------------------------------------------------------------------*/

class pimpleFoamRASState
:
    public simpleFoamState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        pimpleFoamRASState(const pimpleFoamRASState&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleFoamRASState&);


public:

    //- Runtime type information
        TypeName("pimpleFoamRAS");

    // Constructors

        //- Construct from components
        pimpleFoamRASState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~pimpleFoamRASState()
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
