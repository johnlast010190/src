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
    Foam::CassonMu

Description
    An incompressible CassonMu non-Newtonian viscosity model.

    References:
    \verbatim
        CassonMu, N. (1959).
        Rheology of disperse systems.
        In Proceedings of a Conference Organized by the
        British Society of Rheology.
        Pergamon Press, New York.

        Fournier, R. L. (2011).
        Basic transport phenomena in biomedical engineering.
        CRC Press.
    \endverbatim

    Example specification for blood:
    \verbatim
         muModel  CassonMu;

         muModelCoeffs
         {
             m       3.934986e-6;
             tau0    2.9032e-6;
             muMax   13.3333e-6;
             muMin   3.9047e-6;
         }
    \endverbatim

SourceFiles
    CassonMu.C

\*---------------------------------------------------------------------------*/

#ifndef CassonMu_H
#define CassonMu_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CassonMu Declaration
\*---------------------------------------------------------------------------*/

class CassonMu
:
    public materialModel,
    public muModel
{
    // Private data

        //- Model coefficients
        scalar m_;
        scalar tau0_;
        scalar muMin_;
        scalar muMax_;

        //- Quick model access
        enum {sRateInd, modelsEnumSize_};


public:

    //- Runtime type information
    TypeName("CassonMu");


    // Constructors

        //- Construct from components
        CassonMu
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~CassonMu()
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
