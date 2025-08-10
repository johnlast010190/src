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
    Foam::multiComponentFluidState

Description
    Abstract, multiComponentFluid state transport/thermo assembly

    This has been replaced by multiphaseMultiComponentFluidState for USF states
    and should be phased out as legacy solvers are removed.

SourceFiles
    multiComponentFluidState.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentFluidState_H
#define multiComponentFluidState_H

#include "stateComponents/singlePhaseFluid/singlePhaseFluidState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiComponentFluidState Declaration
\*---------------------------------------------------------------------------*/

class multiComponentFluidState
:
    public singlePhaseFluidState
{

protected:

    // Protected Data
    word region_;


public:

    // Static member functions

        //- Single phase multi-species compressible thermophyscial properties
        static Xfer<dictionary> mixtureProperties
        (
            const word& regionName,
            const word& inertSpecie,
            const dictionary& mat
        );

        //- Single phase multi-species compressible species properties
        static Xfer<dictionary> thermoSpeciesProperties
        (
            const wordList& specieNames,
            const dictionary& mat
        );

        //- Single phase multi-species compressible material properties
        static Xfer<dictionary> materialSpeciesProperties
        (
            const wordList& specieNames,
            const dictionary& mat
        );

        //- Single phase multi-species compressible reaction properties
        static Xfer<dictionary> rxnProperties
        (
            const wordList& speciesNames,
            const dictionary& rxnInput
        );

        //- Multi-species combustion properties
        static Xfer<dictionary> combustionProperties
        (
            const wordList& speciesNames,
            const dictionary& combustionInput
        );

    //- Runtime type information
    TypeName("multiComponentFluid");

    // Constructors

        //- Construct from components
        multiComponentFluidState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~multiComponentFluidState()
    {}


    // Member Functions

        //- Assemble field definitions
        virtual void initialise();

        //- Initialise material library
        void initialiseMaterials(dictionary& materialProperties);

        //- Initialise thermo/transport library
        void initialiseThermo(dictionary& materialProperties);

        //- Final manipulations
        virtual void finalise();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
