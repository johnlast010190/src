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
    (c) 2016-2023 Esi Ltd.

Class
    Foam::globalState

Description
    Null state for global stateFunction

SourceFiles
    globalState.C

\*---------------------------------------------------------------------------*/

#ifndef globalState_H
#define globalState_H

#include "stateFunction/stateFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class globalState Declaration
\*---------------------------------------------------------------------------*/

class globalState
:
    public stateFunction
{

private:

    // Private Data

        //- pointer to Time
        const Time* time_;

        //- Disallow default bitwise copy construct
        globalState(const globalState&);

        //- Disallow default bitwise assignment
        void operator=(const globalState&);


public:

    //- Runtime type information
        TypeName("global");

    // Constructors

        //- Construct from components
        globalState
        (
            const Time& time,
            const dictionary& input,
            const dictionary& defaults,
            const bool distributed,
            const bool collated,
            const stateIndex& index
        );

        //- Construct meshed region from components
        globalState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~globalState()
    {}


    // Member Functions

        // Access

            // type indicators

            // region (fluid/solid/surface/none)
            inline virtual regionType region() const
            {
                return rtNone;
            }

            // time (transient/steady)
            inline virtual timeType time() const
            {
                return ttNone;
            }

            // turbulence (RANS/LES/none)
            inline virtual turbulenceType turbulence() const
            {
                return tuNone;
            }

            // compressibility (compressible/incompressible/none)
            inline virtual compressibilityType compressibility() const
            {
                return ctNone;
            }

            const Time& runTime() const;


        // Check

        // Edit

            // Mesh and field construction

                //- generate meshes and store old fields
                virtual void createObjects(const Time& time){}

                //- modify mesh
                virtual void modifyMesh(scalar writePause){}

                //- create fields
                virtual void createFields(){}

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
