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
    (c) 2014-2017 OpenFOAM Foundation
    (c) 2021-2023 Esi Ltd.

Class
    Foam::PengRobinsonGasLaw

Description
    PengRobinsonGas gas equation of state.

SourceFiles
    PengRobinsonGasLawI.H
    PengRobinsonGasLaw.C

\*---------------------------------------------------------------------------*/

#ifndef PengRobinsonGasLaw_H
#define PengRobinsonGasLaw_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PengRobinsonGasLaw Declaration
\*---------------------------------------------------------------------------*/

class PengRobinsonGasLaw
:
    public materialModel,
    public rhoModel,
    public HDepartureModel,
    public CpDepartureModel,
    public SDepartureModel,
    public psiModel,
    public ZModel,
    public CpMCvModel
{
    // Private member data

        //- Gas constant
        scalar R_;

        //- Quick model access
        enum {Rind, Zind, Wind, modelsEnumSize_};

        //- Critical Temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compression factor [-]
        scalar Zc_;

        //- Critical Pressure [Pa]
        scalar Pc_;

        //- Acentric factor [-]
        scalar omega_;

        //- Precalculated coefficients
        scalar a_;
        scalar b_;
        scalar kappa_;


public:

    // Constructors

        //- Construct from dictionary and object registry
        PengRobinsonGasLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<PengRobinsonGasLaw> clone() const;

        //- Selector from dictionary
        inline static autoPtr<PengRobinsonGasLaw> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("PengRobinsonGas");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override;

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override;

        //- Return true if rho is constant
        virtual bool isrhoConst() const override
        {
            return p_->isConst() && T_->isConst();
        }

        //- Return true if Cp departure is constant
        virtual bool isHDepartureConst() const override
        {
            return p_->isConst() && T_->isConst();
        }

        //- Return true if Cp departure is constant
        virtual bool isCpDepartureConst() const override
        {
            return p_->isConst() && T_->isConst();
        }

        //- Return true if Cp departure is constant
        virtual bool isSDepartureConst() const override
        {
            return p_->isConst() && T_->isConst();
        }

        //- Return true if Cp departure is constant
        virtual bool ispsiConst() const override
        {
            return
                T_->isConst()
             && isModelConst<ZModel>() && isModelConst<RModel>();
        }

        //- Return true if CpMCv is constant
        virtual bool isCpMCvConst() const override
        {
            return
                p_->isConst() && T_->isConst()
             && isModelConst<ZModel>() && isModelConst<RModel>();
        }

        //- Functions for density [kg/m^3]
        declareScalarFieldFunctions(rho);

        //- Return enthalpy departure [J/kg]
        declareScalarFieldFunctions(HDeparture);

        //- Return Cp departure [J/(kg K)]
        declareScalarFieldFunctions(CpDeparture);

        //- Functions for entropy [J/(kg K)]
        declareScalarFieldFunctions(SDeparture);

        //- Functions for compressibility rho/p [s^2/m^2]
        declareScalarFieldFunctions(psi);

        //- Functions for compression factor []
        defineScalarFieldFunctionsTwo(Z, p_, T_);

        //- Functions for (Cp - Cv) [J/(kg K)]
        declareScalarFieldFunctions(CpMCv);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update dependent model table
        virtual void updateTable(const word&) override;

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

#include "PengRobinsonGasLawI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
