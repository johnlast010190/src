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
    Foam::BokilBewtraMu

Description
    An incompressible Bokil-Bewtra non-Newtonian viscosity model.

    The model is based on Bokil & Bewtra, 1972. The model blends linear
    to the fluid viscosity under a critical solid concentration value Xc:

    mu = mu_0 * 10^(aX) , X >= Xc
    mu = mu_0 * 10^(aX) * (Xc-X)/Xc + mu_l * (1 - (Xc-X)/Xc) , 0 < X < Xc

SourceFiles
    BokilBewtraMu.C

\*---------------------------------------------------------------------------*/

#ifndef BokilBewtraMu_H
#define BokilBewtraMu_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BokilBewtraMu Declaration
\*---------------------------------------------------------------------------*/

class BokilBewtraMu
:
    public materialModel,
    public muModel
{
    // Private data

        // Concentration field X
        word XName_;

        // Mixture viscosity constant in [kg/ms]
        scalar mu0_;

        // Liquid viscosity in [kg/ms]
        scalar mul_;

        // max viscosity in [kg/ms]
        scalar mumax_;

        // Critical concentration in [g/L]
        scalar Xc_;

        // exponent factor in [m^3/kg]
        scalar a_;

        // used to convert X into [g/L]
        scalar fac_;

        //- Solid tracer field
        const volScalarField* X_;


public:

    //- Runtime type information
    TypeName("BokilBewtraMu");


    // Constructors

        //- Construct from dictionary
        BokilBewtraMu
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    //- Destructor
    virtual ~BokilBewtraMu()
    {}


    // Member Functions

        bool ismuConst() const override
        {
            return false;
        }

        //- Dynamic viscosity [kg/ms]
        declareScalarFieldFunctions(mu);

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
