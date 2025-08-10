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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::kappaPrandtl

Description
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    kappaPrandtl.C

\*---------------------------------------------------------------------------*/

#ifndef kappaPrandtl_H
#define kappaPrandtl_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class kappaPrandtl Declaration
\*---------------------------------------------------------------------------*/

class kappaPrandtl
:
    public materialModel,
    public kappaModel
{
    // Private data

        //- Reciprocal Prandtl Number []
        scalar rPr_;

        //- Quick model access
        enum {Cp, muMod, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary
        kappaPrandtl
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<kappaPrandtl> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("kappaPrandtl");

        //- Is it isotropic model?
        virtual bool isotropic() const override
        {
            return true;
        }

        bool iskappaConst() const override
        {
            return isModelConst<muModel>() && isModelConst<CpModel>();
        }

        //- Thermal conductivity [W/mK]
        declareScalarFieldFunctions(kappa);

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

#endif

// ************************************************************************* //
