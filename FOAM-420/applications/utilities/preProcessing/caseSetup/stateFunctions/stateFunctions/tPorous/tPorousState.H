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
    Foam::tPorousState

Description

SourceFiles
    tPorousState.C

\*---------------------------------------------------------------------------*/

#ifndef tPorousState_H
#define tPorousState_H

#include "stateComponents/porousSolid/porousSolidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tPorousState Declaration
\*---------------------------------------------------------------------------*/

class tPorousState
:
    public porousSolidState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        tPorousState(const tPorousState&);

        //- Disallow default bitwise assignment
        void operator=(const tPorousState&);


public:

    //- Runtime type information
        TypeName("tPorous");

    // Constructors

        //- Construct from components
        tPorousState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~tPorousState()
    {}


    // Member Functions

        // Access

            // type indicators

            // region (fluid/solid/surface/none)
            inline virtual regionType region() const override
            {
                return rtSolid;
            }

            // time (transient/steady)
            inline virtual timeType time() const override
            {
                return ttSteady;
            }

            // turbulence (RANS/LES/none)
            inline virtual turbulenceType turbulence() const override
            {
                return tuNone;
            }

            // compressibility (compressible/incompressible/none)
            inline virtual compressibilityType compressibility() const override
            {
                return ctNone;
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
