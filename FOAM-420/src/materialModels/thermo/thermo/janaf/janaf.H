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
    (c) 2021-2023 Esi Ltd.

Class
    Foam::janaf

Description
    JANAF tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    janafI.H
    janaf.C

\*---------------------------------------------------------------------------*/

#ifndef janaf_H
#define janaf_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class janaf Declaration
\*---------------------------------------------------------------------------*/

class janaf
:
    public materialModel,
    public CpThermoModel,
    public CpModel,
    public HaThermoModel,
    public HaModel,
    public HsModel,
    public HcModel,
    public SThermoModel,
    public SModel,
    public limitModel
{
public:

    // Public data

        static const int nCoeffs_ = 7;
        typedef FixedList<scalar, nCoeffs_> coeffArray;

private:

    // Private member data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;

        //- R
        scalar R_;

        bool isConverted_;

        //- Quick model access
        enum {CpDeparture, HDeparture, SDeparture, HaAll, HcAll, R, modelsEnumSize_};


    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar T) const;

public:


    // Constructors

        //- Construct from dictionary
        janaf
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<janaf> clone() const;

        //- Selector from dictionary
        inline static autoPtr<janaf> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("janaf");

        virtual bool isCpConst() const override
        {
            return T_->isConst() && isModelConst<CpDepartureModel>();
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

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        defineScalarFieldFunctionsOne(limit, T_);

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

#include "janafI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
