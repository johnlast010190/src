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
    Foam::multiphaseFluidState

Description
    Abstract, multiphase state common functions

SourceFiles
    multiphaseFluidState.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseFluidState_H
#define multiphaseFluidState_H

#include "stateComponents/turbulenceModel/turbulenceModelState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stateFunctions
{

/*---------------------------------------------------------------------------*\
                          Class multiphaseFluidState Declaration
\*---------------------------------------------------------------------------*/

class multiphaseFluidState
:
    public turbulenceModelState
{

protected:

    // Protected Data

        wordList phases_;

    // Protected member functions

        //- translate generic phase field gradient entry to phase specific entry
        void addAlphaTypeGrad(const word& phaseName);

        //- binary data assembly, normally scalar
        template<class Type>
        Xfer<HashTable<Type, interfacePair, interfacePair::hash>>
        assembleBinaryData
        (
            word property,
            const Type& defaultValue,
            const List<dictionary>&
        );

        //- binary data assembly, dictionaries
        Xfer<HashTable<dictionary, interfacePair, interfacePair::hash>>
        assembleBinaryDictionaries
        (
            word property,
            const dictionary& defaultValue,
            const List<dictionary>&
        );

        //- combine defaults with user input for material properties
        Xfer<List<dictionary>> assembleMaterialProperties();


private:

    // Private Data


    // Private Functions

        //- Disallow default bitwise copy construct
        multiphaseFluidState(const multiphaseFluidState&);

        //- Disallow default bitwise assignment
        void operator=(const multiphaseFluidState&);



public:

    //- Runtime type information
        TypeName("multiphaseFluid");

    // Constructors

        //- Construct from components
        multiphaseFluidState
        (
            word region,
            const dictionary& input,
            const dictionary& defaults,
            const stateFunction& master,
            const stateIndex& index,
            word meshName = word::null
        );


    //- Destructor
    virtual ~multiphaseFluidState()
    {}


    // Member Functions

        // Access

        // Check

        // Edit

            //- remove default alpha entry from fieldMaps
            virtual void correct();

            //- final manipulations
            virtual void finalise();

        // Write

    // Member Operators


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "multiphaseFluidStateTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
