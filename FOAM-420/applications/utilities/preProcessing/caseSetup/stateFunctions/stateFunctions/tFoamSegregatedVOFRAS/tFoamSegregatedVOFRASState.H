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
    Foam::tFoamSegregatedVOFRASState

Description
    interFoam state identifiers and specialised assembly functions

SourceFiles
    tFoamSegregatedVOFRASState.C

\*---------------------------------------------------------------------------*/

#ifndef tFoamSegregatedVOFRASState_H
#define tFoamSegregatedVOFRASState_H

#include "stateComponents/VOF/VOFState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class tFoamSegregatedVOFRASState Declaration
\*---------------------------------------------------------------------------*/

class tFoamSegregatedVOFRASState
:
    public VOFState
{

private:

    // Private Data

        //- Disallow default bitwise copy construct
        tFoamSegregatedVOFRASState(const tFoamSegregatedVOFRASState&);

        //- Disallow default bitwise assignment
        void operator=(const tFoamSegregatedVOFRASState&);


public:

    //- Runtime type information
        TypeName("tFoamSegregatedVOFRAS");

    // Constructors

        //- Construct from components
        tFoamSegregatedVOFRASState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~tFoamSegregatedVOFRASState()
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
                return tuRAS;
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
