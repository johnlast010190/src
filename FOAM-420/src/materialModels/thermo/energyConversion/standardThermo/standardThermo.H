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
    Foam::standardThermo

SourceFiles
    standardThermo.C

\*---------------------------------------------------------------------------*/

#ifndef standardThermo_H
#define standardThermo_H

#include "materialModels/materialTables/materialTables.H"
#include "materialModels/materialModel/materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class standardThermo Declaration
\*---------------------------------------------------------------------------*/

class standardThermo
:
    public materialModel,
    public AModel,
    public CvModel,
    public EaModel,
    public EsModel,
    public gammaModel,
    public GModel
{
    // Private data

        //- Quick model access
        enum {rho, CpMCv, Cp, S, Ha, Hs, Ea, modelsEnumSize_};


public:

    //- Runtime type information
    TypeName("standardThermo");


    // Constructors

        //- Construct from components
        standardThermo
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- New constructor
        static autoPtr<standardThermo> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    ~standardThermo();

    //- Model dependence
    static const List<List<bool>> dependencies_;


     // Member functions

        virtual bool isAConst() const override
        {
            return
                T_->isConst()
             && isModelConst<EaModel>() && isModelConst<SModel>();
        }

        virtual bool isCvConst() const override
        {
            return isModelConst<CpModel>() && isModelConst<CpMCvModel>();
        }

        virtual bool isEaConst() const override
        {
            return
                p_->isConst()
             && isModelConst<HaModel>() && isModelConst<rhoModel>();
        }

        virtual bool isEsConst() const override
        {
            return
                p_->isConst()
             && isModelConst<HsModel>() && isModelConst<rhoModel>();
        }

        virtual bool isgammaConst() const override
        {
            return isModelConst<CpModel>() && isModelConst<CpMCvModel>();
        }

        virtual bool isGConst() const override
        {
            return
                T_->isConst()
             && isModelConst<HaModel>() && isModelConst<SModel>();
        }

        //- Heat capacity at constant volume [J/(kg K)]
        declareScalarFieldFunctions(Cv);

        //- Gamma = Cp/Cv []
        declareScalarFieldFunctions(gamma);

        //- Sensible internal energy [J/kg]
        declareScalarFieldFunctions(Es);

        //- Absolute internal energy [J/kg]
        declareScalarFieldFunctions(Ea);

        //- Gibbs free energy [J/kg]
        declareScalarFieldFunctions(G);

        //- Helmholtz free energy [J/kg]
        declareScalarFieldFunctions(A);

        //- Return model dependence list
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

#endif

// ************************************************************************* //
