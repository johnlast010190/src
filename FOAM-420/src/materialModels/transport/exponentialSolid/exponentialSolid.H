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
    (c) 2022-2023 Esi Ltd.

Class
    Foam::exponentialSolid

Description
    Exponential properties for solid heat transport.

SourceFiles
    exponentialSolid.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialSolid_H
#define exponentialSolid_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class exponentialSolid Declaration
\*---------------------------------------------------------------------------*/

class exponentialSolid
:
    public materialModel,
    public kappaModel,
    public vKappaModel
{
    // Private data

        //- Constant thermal coefficient.
        scalar kappa0_;

        //- Exponent coefficient
        scalar n0_;

        //- Reference temperature
        scalar Tref_;


public:

    // Constructors

        //- Construct from dictionary
        exponentialSolid
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<exponentialSolid> clone() const;


    // Member functions

        //- Return the instantiated type name
        TypeName("exponentialSolid");

        //- Is it isotropic model?
        virtual bool isotropic() const override
        {
            return true;
        }

        bool iskappaConst() const override
        {
            return T_->isConst();
        }

        bool isvKappaConst() const override
        {
            return T_->isConst();
        }

        //- Isotropic thermal conductivity [W/mK]
        declareScalarFieldFunctions(kappa);

        //- Un-isotropic thermal conductivity [W/mK]
        declareVectorFieldFunctions(vKappa);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update table
        virtual void updateTable(const word& modelName) override
        {}

        //- Models casting
        virtual baseModels<scalar>* castScalarModel
        (
            const word& modelName
        ) override;

        virtual baseModels<vector>* castVectorModel
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
