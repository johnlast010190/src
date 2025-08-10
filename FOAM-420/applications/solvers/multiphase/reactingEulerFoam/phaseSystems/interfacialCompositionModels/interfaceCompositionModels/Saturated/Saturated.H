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
    (c) 2015-2016 OpenFOAM Foundation

Class
    Foam::Saturated

Description
    Model which uses a saturation pressure model for a single species to
    calculate the interface composition.

SourceFiles
    Saturated.C

\*---------------------------------------------------------------------------*/

#ifndef Saturated_H
#define Saturated_H

#include "interfacialCompositionModels/interfaceCompositionModels/InterfaceCompositionModel/InterfaceCompositionModel.H"
#include "interfacialCompositionModels/saturationModels/saturationModel/saturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                         Class Saturated Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class Saturated
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
protected:

    // Private data

        //- Saturated species name
        word saturatedName_;

        //- Saturated species index
        label saturatedIndex_;

        //- Saturation pressure model
        autoPtr<saturationModel> saturationModel_;


    // Private Member Functions

        //- Constant of propotionality between partial pressure and mass
        //  fraction
        tmp<volScalarField> wRatioByP() const;


public:

    //- Runtime type information
    TypeName("Saturated");

    // Constructors

        //- Construct from components
        Saturated
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Saturated();


    // Member Functions

        //- Update the composition
        virtual void update(const volScalarField& Tf);

        //- The interface species fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- The interface species fraction derivative w.r.t. temperature
        virtual tmp<volScalarField> YfPrime
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceCompositionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Saturated.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
