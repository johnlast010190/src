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
    Foam::singlePhaseFluidState

Description
    Abstract, singlePhaseFluid state transport/thermo assembly

SourceFiles
    singlePhaseFluidState.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseFluidState_H
#define singlePhaseFluidState_H

#include "stateComponents/turbulenceModel/turbulenceModelState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class singlePhaseFluidState Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseFluidState
:
    public turbulenceModelState
{


protected:
    // Protected Data


    // Protected member functions


private:

    // Private Data

        word materialName_;

        wordList componentNames_;

    // Private Functions

        //- Disallow default bitwise copy construct
        singlePhaseFluidState(const singlePhaseFluidState&);

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseFluidState&);



public:

    //- static member functions

        //- single phase incompressible transport properties
        static Xfer<dictionary> incompressibleProperties
        (
            const word& materialName,
            const dictionary& materialProperties
        );

        //- single phase compressible thermophyscial properties
        static Xfer<dictionary> compressibleProperties
        (
            const word& materialName,
            const dictionary& materialProperties
        );


    //- Runtime type information
        TypeName("singlePhaseFluid");

    // Constructors

        //- Construct from components
        singlePhaseFluidState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~singlePhaseFluidState()
    {}


    // Member Functions

        // Access

            //const reference to name of fluid
            inline const word& materialName() const
            {
                return materialName_;
            }

            //const reference to name of fluid
            inline const wordList& componentNames() const
            {
                return componentNames_;
            }

        // Check

        // Edit

            //- assemble field definitions
            virtual void initialise();

            //- single phase compressible thermophyscial properties
            Xfer<dictionary> setupMaterialProperties
            (
                const word& materialName,
                const dictionary& materialProperties
            );

            //- Initialise material library
            void initialiseMaterials
            (
                dictionary& materialProperties,
                const dictionary& matDefaults,
                const dictionary& matInputs,
                const bool isMatInInput
            );

            //- Initialise thermo library
            void initialiseThermo
            (
                dictionary& materialProperties,
                const dictionary& matDefaults,
                const dictionary& matInputs,
                const bool isMatInInput
            );

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
