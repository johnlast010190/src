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
    Foam::turbulenceModelState

Description
    turbulenceModel state identifiers and specialised assembly functions

SourceFiles
    turbulenceModelState.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceModelState_H
#define turbulenceModelState_H

#include "stateComponents/region/regionState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class turbulenceModelState Declaration
\*---------------------------------------------------------------------------*/

class turbulenceModelState
:
    public regionState
{


protected:

    // Protected Data

        wordList turbulenceFields_;

    // Protected member functions

        virtual primitiveEntry divScheme(const word&, const dictionary&) const;

        static void addTurbulenceSpecificSchemes(const dictionary& turbDict, dictionary& fvSchemes);


private:

    // Private Data

        //- Disallow default bitwise copy construct
        turbulenceModelState(const turbulenceModelState&);

        //- Disallow default bitwise assignment
        void operator=(const turbulenceModelState&);


    // Private Functions

        //- add solver and relaxation entries
        void fvSolutionInject
        (
            const dictionary& fieldDict,
            const dictionary& modelDict,
            dictionary& fvSolDict,
            const word& tfname
        );


public:

    //- Runtime type information
        TypeName("turbulenceModel");

    // Constructors

        //- Construct from components
        turbulenceModelState
        (
            const word& region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            const word& meshName = word::null
        );


    //- Destructor
    ~turbulenceModelState() override = default;


    // Member Functions

        // Access

            // type indicators

            // region (fluid/solid/surface/none)
            inline regionType region() const override
            {
                return rtFluid;
            }

        // Check

        // Edit

            //- assemble field definitions
            void initialise() override;

            //- correct stateDict
            void correct() override;

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
