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
    (c) 2016 OpenFOAM Foundation
    (c) 2021-2023 Esi Ltd.

Class
    Foam::strainRateFunctionMu

Description
    Run-time selected strain-rate function non-Newtonian viscosity model.

    Example linear function of strain-rate:
    \verbatim
        muModel  strainRateFunctionMu;

        muModelCoeffs
        {
            function polynomial ((0 0.1) (1 1.3));
        }
    \endverbatim

See also
    Foam::Function1

SourceFiles
    strainRateFunctionMu.C

\*---------------------------------------------------------------------------*/

#ifndef strainRateFunctionMu_H
#define strainRateFunctionMu_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"
#include "primitives/functions/Function1/Function1/Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class strainRateFunctionMu Declaration
\*---------------------------------------------------------------------------*/

class strainRateFunctionMu
:
    public materialModel,
    public muModel
{
    // Private data

        //- Strain-rate function
        autoPtr<Function1<scalar>> strainRateFunctionMu_;

        //- Quick model access
        enum {sRateInd, modelsEnumSize_};


public:

    //- Runtime type information
    TypeName("strainRateFunctionMu");


    // Constructors

        //- Construct from components
        strainRateFunctionMu
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~strainRateFunctionMu()
    {}


    // Member Functions

        bool ismuConst() const override
        {
            return isModelConst<strainRateModel>();
        }

        //- Dynamic viscosity [kg/ms]
        declareScalarFieldFunctions(mu);

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

#endif

// ************************************************************************* //
