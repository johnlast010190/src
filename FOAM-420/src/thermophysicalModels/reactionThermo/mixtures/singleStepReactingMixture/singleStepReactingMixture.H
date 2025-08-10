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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::singleStepReactingMixture

Group
    grpReactionThermophysicalMixtures

Description
    Single step reacting mixture

SourceFiles
    singleStepReactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef singleStepReactingMixture_H
#define singleStepReactingMixture_H

#include "chemistryReaders/chemistryReader/chemistryReader.H"
#include "mixtures/reactingMixture/reactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class singleStepReactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class singleStepReactingMixture
:
    public reactingMixture<ThermoType>
{

protected:

    // Protected data

        //- Stoichiometric air-fuel mass ratio
        dimensionedScalar stoicRatio_;

        //- Stoichiometric oxygen-fuel mass ratio
        dimensionedScalar s_;

        //- Heat of combustion [J/Kg]
        dimensionedScalar qFuel_;

        //- Stoichiometric coefficient for the reaction.
        scalarList specieStoichCoeffs_;

        //- Mass concentrations at stoichiometric mixture for fres.
        scalarList Yprod0_;

        //- List of components residual
        PtrList<volScalarField> fres_;

        //- Inert specie index
        label inertIndex_;

        //- Fuel specie index
        label fuelIndex_;

        //- List to indicate if specie is produced/consumed
        List<int> specieProd_;


    // Protected member functions

        //- Calculate qFuel
        void calculateqFuel();

        //- Calculate air/fuel and oxygen/fuel ratio
        void massAndAirStoichRatios();

        //- Calculate maximum products at stoichiometric mixture
        void calculateMaxProducts();


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        singleStepReactingMixture(const singleStepReactingMixture&);

        //- Disallow default bitwise assignment
        void operator=(const singleStepReactingMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        singleStepReactingMixture
        (
            const dictionary&,
            const objectRegistry&,
            const word&
        );


    //- Destructor
    virtual ~singleStepReactingMixture()
    {}


    // Member functions

        //- Calculates the residual for all components
        void fresCorrect();


        // Access functions

            //- Return the stoichiometric air-fuel mass ratio
            inline const dimensionedScalar stoicRatio() const;

            //- Return the Stoichiometric oxygen-fuel mass ratio
            inline const dimensionedScalar s() const;

            //- Return the heat of combustion [J/Kg]
            inline const dimensionedScalar qFuel() const;

            //- Return the stoichiometric coefficient for the reaction
            inline const List<scalar>& specieStoichCoeffs() const;

            //- Return the list of components residual
            inline tmp<volScalarField> fres(const label index) const;

            //- Return the inert specie index
            inline label inertIndex() const;

            //- Return the fuel specie index
            inline label fuelIndex() const;

            //- Return the list to indicate if specie is produced/consumed
            inline const List<int>& specieProd() const;

            //- Return the list of products mass concentrations
            inline const scalarList& Yprod0() const;


        // I-O

            //- Read dictionary
            void read(const objectRegistry&, const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixtures/singleStepReactingMixture/singleStepReactingMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtures/singleStepReactingMixture/singleStepReactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
