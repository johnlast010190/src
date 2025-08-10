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
    Foam::adjointMultipointState

Description
    adjointMultipoint state identifiers and specialised assembly functions

SourceFiles
    adjointMultipointState.C

\*---------------------------------------------------------------------------*/

#ifndef adjointMultipointState_H
#define adjointMultipointState_H

#include "stateFunctions/rhoSimpleFoam/rhoSimpleFoamState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class adjointMultipointState Declaration
\*---------------------------------------------------------------------------*/

class adjointMultipointState
:
    public rhoSimpleFoamState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        adjointMultipointState(const adjointMultipointState&);

        //- Disallow default bitwise assignment
        void operator=(const adjointMultipointState&);



public:

    //- Runtime type information
        TypeName("adjointMultipoint");

    // Constructors

        //- Construct from components
        adjointMultipointState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~adjointMultipointState()
    {}


    // Member Functions

        // Access

            // type indicators

            // region (fluid/solid/surface/none)
            inline regionType region() const override
            {
                return rtMultipointAdjoint;
            }

            //- USF vs legacy solver state
            virtual bool usf() const override
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
