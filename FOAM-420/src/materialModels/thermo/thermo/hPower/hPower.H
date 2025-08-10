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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021-2023 Esi Ltd.

Description
    Power-function based thermodynamics package templated on EquationOfState.

    In this thermodynamics package the heat capacity is a simple power of
    temperature:

        Cp(T) = c0*(T/Tref)^n0;

    which is particularly suitable for solids.

Class
    Foam::hPower

SourceFiles
    hPowerI.H
    hPower.C

\*---------------------------------------------------------------------------*/

#ifndef hPower_H
#define hPower_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hPower Declaration
\*---------------------------------------------------------------------------*/

class hPower
:
    public materialModel,
    public CpThermoModel,
    public CpModel,
    public HaThermoModel,
    public HaModel,
    public HsThermoModel,
    public HsModel,
    public HcModel,
    public SThermoModel,
    public SModel
{
    // Private member data

        scalar c0_;
        scalar n0_;
        scalar Tref_;
        scalar Hf_;

        //- Quick model access
        enum {HsInd, HcInd, CpDeparture, HDeparture, SDeparture, modelsEnumSize_};


public:


    // Constructors

        //- Construct from dictionary
        hPower
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<hPower> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hPower> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hPower");

        virtual bool isCpConst() const override
        {
            return T_->isConst() && isModelConst<CpDepartureModel>();
        }

        virtual bool isHaThermoConst() const override
        {
            return isModelConst<HsModel>() && isModelConst<HcModel>();
        }

        virtual bool isHaConst() const override
        {
            return isHaThermoConst() && isModelConst<HDepartureModel>();
        }

        virtual bool isHsConst() const override
        {
            return T_->isConst() && isModelConst<HDepartureModel>();
        }

        virtual bool isSConst() const override
        {
            return T_->isConst() && isModelConst<SDepartureModel>();
        }

        //- Heat capacity at constant pressure [J/(kg K)]
        defineScalarFieldFunctionsOne(CpThermo, T_);

        //- Absolute Enthalpy [J/kg]
        declareScalarFieldFunctions(HaThermo);

        //- Sensible enthalpy [J/kg]
        defineScalarFieldFunctionsOne(HsThermo, T_);

        //- Entropy [J/(kg K)]
        defineScalarFieldFunctionsOne(SThermo, T_);

        //- Heat capacity at constant pressure [J/(kg K)]
        declareScalarFieldFunctions(Cp);

        //- Absolute Enthalpy [J/kg]
        declareScalarFieldFunctions(Ha);

        //- Sensible enthalpy [J/kg]
        declareScalarFieldFunctions(Hs);

        //- Chemical enthalpy [J/kg]
        defineScalarFieldFunctionsZero(Hc);

        //- Entropy [J/(kg K)]
        declareScalarFieldFunctions(S);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update table
        virtual void updateTable(const word& modelName) override;

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

#include "hPowerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
