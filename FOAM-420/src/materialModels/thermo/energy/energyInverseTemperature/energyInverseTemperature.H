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
    (c) 2021 Esi Ltd.

Class
    energyInverseTemperature

SourceFiles
    energyInverseTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef energyInverseTemperature_H
#define energyInverseTemperature_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/materialTables/materialTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class energyInverseTemperature Declaration
\*---------------------------------------------------------------------------*/

class energyInverseTemperature
:
    public materialModel,
    public TModel,
    public THaModel,
    public THsModel,
    public TEaModel,
    public TEsModel
{

    // Static data members

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;


    // Private data members

        //- Temperature
        volScalarField* TnonConst_;

        //- Reference offset
        scalar Tref_;

        //- Energy
        const volScalarField* he_;

        //- Quick model access
        enum {HE, Cpv, limit, Ha, Hs, Ea, Es, Cp, Cv, modelsEnumSize_};


    // Private member functions

        //- Lookup energy field name
        word lookupEnergyFieldName();


public:


    // Constructors

        //- Construct from components
        energyInverseTemperature
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~energyInverseTemperature();


    // Member functions

        //- Return the instantiated type name
        TypeName("energyInverse");

        //- Mixture temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(p, T)
        //  and dF(p, T)/dT
        // For cells
        scalar Cell
        (
            const label celli,
            const label HEFun,
            const label CpvFun
        ) const;

        //- For patches
        tmp<scalarField> Patch
        (
            const label patchi,
            const label HEFun,
            const label CpvFun
        ) const;

        //- For internal fields
        tmp<scalarField> Internal
        (
            const label HEFun,
            const label CpvFun
        ) const;

        bool isTConst() const override
        {
            return T_->isConst();
        }

        bool isTHsConst() const override
        {
            return T_->isConst();
        }

        bool isTHaConst() const override
        {
            return T_->isConst();
        }

        bool isTEsConst() const override
        {
            return T_->isConst();
        }

        bool isTEaConst() const override
        {
            return T_->isConst();
        }

        //- Inversion based on selected energy model
        declareScalarFieldFunctions(T);

        //- Inversion with Ha
        declareScalarFieldFunctions(THa);

        //- Inversion with Hs
        declareScalarFieldFunctions(THs);

        //- Inversion with Ea
        declareScalarFieldFunctions(TEa);

        //- Inversion with Es
        declareScalarFieldFunctions(TEs);

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

        //- Update fields
        virtual bool updateScalarField
        (
            const word& fieldName,
            const volScalarField& volField
        ) override;

        //- Re-read data
        virtual bool read() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
