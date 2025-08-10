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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2021-2023 Esi Ltd.

Class
    Foam::matHeThermo

Description
    Enthalpy/Internal energy for a mixture

SourceFiles
    matHeThermo.C

\*---------------------------------------------------------------------------*/

#ifndef matHeThermo_H
#define matHeThermo_H

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
                         Class matHeThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class MixtureType>
class matHeThermo
:
    public BasicThermo,
    public MixtureType
{
protected:

    // Protected data

        //- Energy
        volScalarField he_;

        //- Reference to material models
        materialTables& materials_;

        //- Is the fluid buoyant?
        bool isBuoyant_;

        //- Is boyancy distinct?
        bool isDistinctBuoyancy_;

        //- Do we solve for T from he or vice versa?
        bool solveTFromhe_;


    // Protected Member Functions

        // Enthalpy/Internal energy

            //- Correct the enthalpy/internal energy field boundary
            void heBoundaryCorrection(fvPatchScalarField& pf);

            //- Correct the enthalpy/internal energy field boundaries
            void heBoundaryCorrection(volScalarField& he);

            //- Construct material tables
            materialTables& matLookupOrConstruct
            (
                const objectRegistry& obr,
                const dictionary& dict
            );


private:

    // Private Member Functions

        //- Initialize
        void init();

        //- Initialize all the he times
        void initHe
        (
            const volScalarField& p,
            volScalarField& he
        );

        //- Update T, p
        tmp<scalarField> tmpBoundaryField
        (
            const scalarField& p,
            const scalarField& T,
            const label patchi,
            const word& modelName
        ) const;

        //- Create name of the energy field based on the model
        word energyName(const dictionary& dict) const;

        //- Initialize buoyancy flags
        void initBuoyancyFlags(const word& buoyantName);

        //- Report solution flags
        void reportFlags() const;


public:

    // Constructors

        //- Construct from object registry, dictionary and phase name
        matHeThermo
        (
            const objectRegistry&,
            const dictionary&,
            const word& phaseName = word::null
        );

        //- Construct from object registry and phase name
        matHeThermo
        (
            const objectRegistry&,
            const word& phaseName = word::null
        );

        //- Return material tables
        virtual materialTables& materials() const
        {
            return materials_;
        }




    //- Destructor
    virtual ~matHeThermo();


    // Member functions

        //- Return the compostion of the mixture
        virtual typename MixtureType::basicMixtureType& composition()
        {
            return *this;
        }

        //- Return the compostion of the mixture
        const virtual typename MixtureType::basicMixtureType&
        composition() const
        {
            return *this;
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const;

        //- Return true if Cp or Cv (depending selected energy type) is constant
        virtual bool isCpvConst() const;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const;

        //- Is the model buoyant?
        virtual bool buoyant() const
        {
            return isBuoyant_;
        }

        //- It the buoyant rho same as main rho?
        virtual bool distinctBuoyancy() const
        {
            return isDistinctBuoyancy_;
        }

        //- Do we solve for T from he or vice-versa?
        virtual bool solvesTFromhe() const
        {
            return solveTFromhe_;
        }

        //- Set to solve for T from he
        virtual void setSolveTfromhe()
        {
            solveTFromhe_ = true;
        }

        //- Set to solve for he from T
        virtual void setSolveheFromT()
        {
            solveTFromhe_ = false;
        }


        // Access to thermodynamic state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                return he_;
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                return he_;
            }


        // Fields derived from thermodynamic state variables

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& p,
                const scalarField& T0,
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& p,
                const scalarField& T0,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const;


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const;

            //- Thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const;

            //- Effective thermal diffusivity for temperature
            //  of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff(const volScalarField&)
            const;

            //- Effective thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Density [kg/m^3] (density for use in buoyancy term)
            virtual tmp<volScalarField> buoyantRho() const;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "matHeThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
