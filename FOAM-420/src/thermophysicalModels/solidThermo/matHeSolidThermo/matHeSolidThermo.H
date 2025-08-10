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
    (c) 2022 Esi Ltd.
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::matHeSolidThermo

Description
    Energy for a solid mixture

SourceFiles
    matHeSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef matHeSolidThermo_H
#define matHeSolidThermo_H

#include "materials/matHeThermo/matHeThermo.H"
#include "materialModels/materialTables/materialTables.H"
#include "materialsSpecie/specieModels/specieModels.H"
#include "equationOfState/equationOfState/equationOfState.H"
#include "thermo/thermo/thermodynamics/thermodynamics.H"
#include "thermo/energyConversion/energyConversion/energyConversion.H"
#include "thermo/energy/energy/energy.H"
#include "thermo/energy/energyInverseTemperature/energyInverseTemperature.H"
#include "thermo/energyConversion/standardThermo/standardThermo.H"
#include "materialsSpecie/materialsSpecie.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class matHeSolidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicSolidThermo, class MixtureType>
class matHeSolidThermo
:
    public matHeThermo<BasicSolidThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        matHeSolidThermo(const matHeSolidThermo<BasicSolidThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("matHeSolidThermo");


    // Constructors

        //- Construct from obr and phase name
        matHeSolidThermo
        (
            const objectRegistry&,
            const word& phaseName
        );

        //- Construct from obr, dictionary and phase name
        matHeSolidThermo
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName
        );


    //- Destructor
    virtual ~matHeSolidThermo();


    // Member functions

        //- Update properties
        virtual void correct();


        // Derived thermal properties

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const;

            //- Return true if thermal conductivity is isotropic
            virtual bool isotropic() const
            {
                return matHeThermo<BasicSolidThermo, MixtureType>::isotropic();
            }


        // Per patch calculation

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<vectorField> Kappa(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "matHeSolidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
