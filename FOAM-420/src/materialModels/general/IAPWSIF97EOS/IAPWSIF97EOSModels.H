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

Class
    Foam::IAPWSIF97EOSModels

Description
    This calculation page is based on the Thermodynamic Properties of Water
    and Steam (Reference) [1] (named as IAPWS-IF97) provided by IAPWS.
    Detailed information about used equations, constants, range of validity
    etc is presented in PDF version of IAPWS.

    Release which can be downloaded from IAPWS web site www.iapws.org

Note:
    1. These models are not fully consistent with the material library.
    Additional work is required to fully integrate.
    Now the gas constant R is hard-coded + in case of computing Cp
    the rho model is requred for region3 that is used directly from
    standard.
    2. No departure models used.

Reference
    1. Revised Release on the IAPWS Industrial Formulation 1997 for
    the Thermodynamic Properties of Water and Steam.
    Lucerne, Switzerland, August 2007, available at http://www.iapws.org
    2. Kretzschmar H, Wagner W International Steam Tables.

SourceFiles
    IAPWSIF97EOSModelsI.H
    IAPWSIF97EOSModels.C

\*---------------------------------------------------------------------------*/

#ifndef IAPWSIF97EOSModels_H
#define IAPWSIF97EOSModels_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"
#include "general/IAPWSIF97EOS/IAPWSIF97EOS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IAPWSIF97EOSModels Declaration
\*---------------------------------------------------------------------------*/

class IAPWSIF97EOSModels
:
    public materialModel,
    public rhoModel,
    public psiModel,
    public ZModel,
    public CpMCvModel,
    public CpModel,
    public CvModel,
    public HaModel,
    public HsModel,
    public HcModel,
    public SModel,
    public muModel,
    public kappaModel
{
    // Private member data

        //- Model
        IAPWSIF97EOS waterSteamModel_;

        //- Chemical enthalpy [J/kg]
        scalar Hf_;

        //- Entropy [J/(kg K)]
        scalar Sf_;

        //- Should the model be allowed to cross to saturation curve?
        bool supportPhaseChange_;

        //- Region number
        label regionNumber_;

        //- Checking if we are not crossing to other regions
        inline void checkPhaseChange
        (
            const scalar p,
            const scalar T
        ) const;


public:

    // Constructors

        //- Construct from dictionary and object registry
        IAPWSIF97EOSModels
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<IAPWSIF97EOSModels> clone() const;

        //- Selector from dictionary
        inline static autoPtr<IAPWSIF97EOSModels> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("IAPWSIF97EOS");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override
        {
            return false;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override
        {
            return false;
        }

        //- Functions for density [kg/m^3]
        defineScalarFieldFunctionsTwo(rho, p_, T_);

        //- Functions for compressibility rho/p [s^2/m^2]
        defineScalarFieldFunctionsTwo(psi, p_, T_);

        //- Functions for compression factor []
        defineScalarFieldFunctionsNone(Z);

        //- Functions for (Cp - Cv) [J/(kg K)]
        defineScalarFieldFunctionsTwo(CpMCv, p_, T_);

        //- Absolute Enthalpy [J/kg]
        defineScalarFieldFunctionsTwo(Ha, p_, T_);

        //- Sensible enthalpy [J/kg]
        defineScalarFieldFunctionsTwo(Hs, p_, T_);

        //- Chemical enthalpy [J/kg]
        defineScalarFieldFunctionsZero(Hc);

        //- Entropy [J/(kg K)]
        defineScalarFieldFunctionsTwo(S, p_, T_);

        //- Functions for (Cp - Cv) [J/(kg K)]
        defineScalarFieldFunctionsTwo(Cp, p_, T_);

        //- Functions for Cv [J/(kg K)]
        defineScalarFieldFunctionsTwo(Cv, p_, T_);

        //- Dynamic viscosity [kg/ms]
        defineScalarFieldFunctionsTwo(mu, p_, T_);

        //- Thermal conductivity [W/mK]
        defineScalarFieldFunctionsTwo(kappa, p_, T_);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update dependent model table
        virtual void updateTable(const word&) override
        {}

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IAPWSIF97EOSModelsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
