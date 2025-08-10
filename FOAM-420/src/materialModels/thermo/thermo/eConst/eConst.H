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
    Foam::eConst

SourceFiles
    eConstI.H
    eConst.C

\*---------------------------------------------------------------------------*/

#ifndef eConst_H
#define eConst_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class eConst Declaration
\*---------------------------------------------------------------------------*/

class eConst
:
    public materialModel,
    public CpModel,
    public HaModel,
    public HsModel,
    public HcModel,
    public SModel
{
    // Private member data

        //- Heat capacity [J/(kg K)]
        dimensionedScalar Cv_;

        //- Chemical enthalpy [J/kg]
        dimensionedScalar Hf_;

        //- Quick model access
        enum {CpMCv, Cp, CpDeparture, HDeparture, SDeparture, modelsEnumSize_};


public:


    // Constructors

        //- Construct from dictionary
        eConst
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<eConst> clone() const;

        //- Selector from dictionary
        inline static autoPtr<eConst> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member Functions

        //- Return the instantiated type name
        TypeName("eConst");

        virtual bool isCpConst() const override
        {
            return
                isModelConst<CpMCvModel>() && isModelConst<CpDepartureModel>();
        }

        virtual bool isHaConst() const override
        {
            return
                T_->isConst()
             && isModelConst<CpModel>() && isModelConst<HDepartureModel>();
        }

        virtual bool isHsConst() const override
        {
            return
                T_->isConst()
             && isModelConst<CpModel>() && isModelConst<HDepartureModel>();
        }

        virtual bool isSConst() const override
        {
            return T_->isConst() && isModelConst<SDepartureModel>();
        }

        //- Chemical enthalpy [J/kg]
        defineScalarFieldFunctionsZero(Hc);

        //- Heat capacity at constant pressure [J/(kg K)]
        declareScalarFieldFunctions(Cp);

        //- Absolute Enthalpy [J/kg]
        declareScalarFieldFunctions(Ha);

        //- Sensible enthalpy [J/kg]
        declareScalarFieldFunctions(Hs);

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

#include "eConstI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
