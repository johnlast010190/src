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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019-2023 Esi Ltd.

Class
    Foam::hFunctionT
Description
    Thermodynamics package templated on the equation of state, using Function1
    functions for \c cp, \c h and \c s.

    Functions for \c h and \c s derived from \c cp.

Usage

    \table
        Property     | Description
        Hf           | Heat of formation
        Sf           | Standard entropy
        Cp           | Specific heat at constant pressure function
    \endtable

    Example of the specification of the thermodynamic properties:
    \verbatim
    thermodynamics
    {
        Hf              0;
        Sf              0;
        Cp              polynomial ( ( 1000 0 ) (-0.05 1) (0.003 2) );
    }
    \endverbatim

    The polynomial expression is evaluated as so:

        \f[
            Cp = 1000 - 0.05 T + 0.003 T^2
        \f]

Note
    - Heat of formation is inputted in [J/kg], but internally uses [J/kmol]
    - Standard entropy is inputted in [J/kg/K], but internally uses [J/kmol/K]
    - Specific heat at constant pressure polynomial coefficients evaluate to an
      expression in [J/(kg.K)].

SourceFiles
    hFunctionTI.H
    hFunctionT.C

\*---------------------------------------------------------------------------*/

#ifndef hFunctionT_H
#define hFunctionT_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"
#include "primitives/functions/Function1/Constant/Constant.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hFunctionT Declaration
\*---------------------------------------------------------------------------*/

class hFunctionT
:
    public materialModel,
    public CpThermoModel,
    public CpModel,
    public HaThermoModel,
    public HaModel,
    public HsModel,
    public HcModel,
    public SThermoModel,
    public SModel
{
    // Private member data

        //- Heat of formation
        scalar Hf_;

        //- Standard entropy
        scalar Sf_;

        //- Specific heat at constant pressure
        autoPtr<Function1<scalar>> Cp_;

        //- Quick model access
        enum {Ha, HcInd, CpDeparture, HDeparture, SDeparture, modelsEnumSize_};


public:


    // Constructors

        //- Construct from dictionary
        hFunctionT
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<hFunctionT> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hFunctionT> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hFunctionT");

        virtual bool isCpConst() const override
        {
            return
                (
                    Cp_().type() == Function1Types::Constant<scalar>::typeName
                 || T_->isConst()
                )
             && isModelConst<CpDepartureModel>();
        }

        virtual bool isHaConst() const override
        {
            return T_->isConst() && isModelConst<HDepartureModel>();
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
        defineScalarFieldFunctionsOne(HaThermo, T_);

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

#include "hFunctionTI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
