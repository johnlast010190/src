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
    Foam::multiphaseMultiComponentState

Description
    Base state for new material library which can handle one or more phases
    each with one or more components

SourceFiles
    multiphaseMultiComponentState.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseMultiComponentState_H
#define multiphaseMultiComponentState_H

#include "stateComponents/turbulenceModel/turbulenceModelState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseMultiComponentState Declaration
\*---------------------------------------------------------------------------*/

class multiphaseMultiComponentState
:
    public turbulenceModelState
{

private:

    // Private Data

        //- The material names as specified in caseSetupDict. Either the phase
        //  materials (for multiphase) or species materials (for single phase)
        //  depending on the derived state
        wordList materialNames_;


    // Private Functions

        //- Disallow default bitwise copy construct
        multiphaseMultiComponentState(const multiphaseMultiComponentState&);

        //- Disallow default bitwise assignment
        void operator=(const multiphaseMultiComponentState&);


protected:
    // Protected Data


    // Protected member functions

        //- Initialise material properties
        void initialiseMaterials();

        //- Initialise phase sub-dict of material properties
        void initialisePhase
        (
            const word& phaseName,
            dictionary& materialProperties,
            const dictionary& inDict
        );

        //- Initialise component sub-dict of material properties
        void initialiseComponent
        (
            const word& componentName,
            dictionary& materialProperties,
            const dictionary& inDict,
            const dictionary& inDictParent
        );

        //- translate generic phase field gradient entry to phase specific entry
        void addAlphaTypeGrad(const word& phaseName);

public:

    // Static member functions

        //- Species properties
        static Xfer<dictionary> speciesProperties
        (
            const word& matType,
            const dictionary& matSpecies
        );

public:

    //- Runtime type information
        TypeName("multiphaseMultiComponent");

    // Constructors

        //- Construct from components
        multiphaseMultiComponentState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~multiphaseMultiComponentState()
    {}


    // Member Functions

        // Access

            //- The list of materials as specified in caseSetupDict - these
            //  are the top-level materials, either species in a single-phase
            //  state or phase materials in a multi-phase state
            virtual wordList materialNames() const
            {
                return wordList(materialNames_);
            }

            //- Override to provide the names of the phases
            virtual wordList phaseNames() const
            {
                // Default is the single-phase case unless
                // 'materialType multiphase' is specified in the input
                if
                (
                    input().subOrEmptyDict
                    (
                        "materialProperties"
                    ).lookupOrDefault
                    (
                        "materialType", word::null
                    ) == "multiphase"
                )
                {
                    return materialNames();
                }
                else
                {
                    return wordList(1, word::null);
                }
            }


        // Check

        // Edit

            //- assemble field definitions
            virtual void initialise();

            //- Final manipulations
            virtual void finalise();

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
