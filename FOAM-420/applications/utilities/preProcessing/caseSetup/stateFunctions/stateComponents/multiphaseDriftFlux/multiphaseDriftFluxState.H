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
    (c) 2021 Esi Ltd.

Class
    Foam::multiphaseDriftFlux

Description
    Abstract, multiphaseDriftFlux state transport/thermo assembly

SourceFiles
    multiphaseDriftFlux.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseDriftFlux_H
#define multiphaseDriftFlux_H

#include "stateComponents/multiphaseFluid/multiphaseFluidState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseDriftFlux Declaration
\*---------------------------------------------------------------------------*/

class multiphaseDriftFluxState
:
    public multiphaseFluidState
{

protected:

    // Protected Data

    // Protected member functions


private:

    // Private Data


    // Private Functions

        //- scrape diameter model for inclusion into phase properties
        Xfer<dictionary> eulerPhaseProperties
        (
            const word& materialName,
            const dictionary&
        ) const;

        //- initialise multiphase, might extend to twophaseEuler in future
        void initialiseMultiPhase();

        //- Disallow default bitwise copy construct
        multiphaseDriftFluxState(const multiphaseDriftFluxState&);

        //- Disallow default bitwise assignment
        void operator=(const multiphaseDriftFluxState&);



public:

    //- Runtime type information
        TypeName("multiphaseDriftFlux");

    // Constructors

        //- Construct from components
        multiphaseDriftFluxState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~multiphaseDriftFluxState()
    {}


    // Member Functions

        // Access

        // Check

        // Edit

            //- assemble field definitions
            virtual void initialise();

            //- remove default U entry from fieldMaps
            virtual void correct();

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
