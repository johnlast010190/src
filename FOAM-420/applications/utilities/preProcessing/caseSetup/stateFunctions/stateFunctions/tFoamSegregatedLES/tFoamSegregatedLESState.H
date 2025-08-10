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
    (c) 2023 Esi Ltd.

Class
    Foam::tFoamSegregatedLESState

Description
    simpleFoam state identifiers and specialised assembly functions

SourceFiles
    tFoamSegregatedLESState.C

\*---------------------------------------------------------------------------*/

#ifndef tFoamSegregatedLESState_H
#define tFoamSegregatedLESState_H

#include "stateComponents/multiphaseMultiComponent/multiphaseMultiComponentState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tFoamSegregatedLESState Declaration
\*---------------------------------------------------------------------------*/

class tFoamSegregatedLESState
:
    public multiphaseMultiComponentState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        tFoamSegregatedLESState(const tFoamSegregatedLESState&);

        //- Disallow default bitwise assignment
        void operator=(const tFoamSegregatedLESState&);



public:

    //- Runtime type information
        TypeName("tFoamSegregatedLES");

    // Constructors

        //- Construct from components
        tFoamSegregatedLESState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~tFoamSegregatedLESState()
    {}


    // Member Functions

        // Access

            // type indicators

            // region (fluid/solid/surface/none)
            inline virtual regionType region() const
            {
                return rtFluid;
            }

            // time (transient/steady)
            inline virtual timeType time() const
            {
                return ttTrans;
            }

            // turbulence (RANS/LES/none)
            inline virtual turbulenceType turbulence() const
            {
                return tuLES;
            }

            // compressibility (compressible/incompressible/none)
            inline virtual compressibilityType compressibility() const
            {
                return ctComp;
            }

            //- USF vs legacy solver state
            virtual bool usf() const
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
