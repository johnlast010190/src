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
    (c) 2021-2023 Esi Ltd.

Class
    Foam::hTabulated

Description
    Tabulated properties thermodynamics package
    templated into the EquationOfState.

    thermodynamics
    {
        //- Heat capacity table as function of P and T (required)
        CpTableCoeffs
        {
            fileName        "$FOAM_CASE/constant/CpTableName.dat";
            outOfBounds     clamp;
        }

        //- Heat of formation (required)
        Hf 1901;
    }

    The Cp tables are described in terms of
    the pressure and temperature
    For example, the table:

                           temperature (K)
                       |  300   450   600
                  -----+-----------------
                  1e5  |   A     B     C
    pressure (Pa) 1.1e5|   D     E     F
                  2e5  |   G     H     I


    Is specified by the following:

       (
            (
                1e5
                (
                    (300    A)
                    (450    B)
                    (600    C)
                )
            )
            (
                1.1e5
                (
                    (300    D)
                    (450    E)
                    (600    F)
                )
            )
            (
                2e5
                (
                    (300    G)
                    (450    H)
                    (600    I)
                )
            )
        )

    All pressures and temperatures must be in ascending order
    in each list.

SourceFiles
    hTabulatedI.H
    hTabulated.C

\*---------------------------------------------------------------------------*/

#ifndef hTabulated_H
#define hTabulated_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
#include "interpolations/interpolation2DTable/interpolation2DTable.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class hTabulated Declaration
\*---------------------------------------------------------------------------*/

class hTabulated
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

        //- 2D look up table Cp = function of pressure and temperature
        //- Specific heat at constant pressure [J/(kg.K)]
        autoPtr<interpolation2DTable<scalar>> CpTable_;

        //- Heat of formation
        scalar Hf_;

        //- Quick model access
        enum {CpDeparture, HDeparture, SDeparture, modelsEnumSize_};


public:


    // Constructors

        //- Construct from dictionary
        hTabulated
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<hTabulated> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hTabulated> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("hTabulated");

        virtual bool isCpConst() const override
        {
            return
                p_->isConst() && T_->isConst()
             && isModelConst<CpDepartureModel>();
        }

        virtual bool isHaConst() const override
        {
            return
                p_->isConst() && T_->isConst()
             && isModelConst<HDepartureModel>();
        }

        virtual bool isHsConst() const override
        {
            return
                p_->isConst() && T_->isConst()
             && isModelConst<HDepartureModel>();
        }

        virtual bool isSConst() const override
        {
            return
                p_->isConst() && T_->isConst()
             && isModelConst<SDepartureModel>();
        }

        //- Heat capacity at constant pressure [J/(kg K)]
        defineScalarFieldFunctionsTwo(CpThermo, p_, T_);

        //- Absolute Enthalpy [J/kg]
        defineScalarFieldFunctionsTwo(HaThermo, p_, T_);

        //- Sensible enthalpy [J/kg]
        defineScalarFieldFunctionsTwo(HsThermo, p_, T_);

        //- Entropy [J/(kg K)]
        defineScalarFieldFunctionsTwo(SThermo, p_, T_);

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

#include "hTabulatedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
