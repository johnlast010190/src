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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2021-2022 Esi Ltd.

Class
    Foam::BoussinesqLaw

Description
    Incompressible gas equation of state using the Boussinesq approximation for
    the density as a function of temperature only:

    \verbatim
        rho = rho0*(1 - beta*(T - T0))
    \endverbatim

SourceFiles
    BoussinesqLawI.H
    BoussinesqLaw.C

\*---------------------------------------------------------------------------*/

#ifndef BoussinesqLaw_H
#define BoussinesqLaw_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BoussinesqLaw Declaration
\*---------------------------------------------------------------------------*/

class BoussinesqLaw
:
    public materialModel,
    public rhoModel,
    public psiModel,
    public ZModel,
    public CpMCvModel
{
    // Private member data

        //- Reference density
        scalar rho0_;

        //- Reference temperature
        scalar T0_;

        //- Thermal expansion coefficient
        scalar beta_;

        //- Gas constant
        scalar R_;

        //- Quick model access
        enum {R, modelsEnumSize_};


public:

    // Constructors

        //- Construct from dictionary and object registry
        BoussinesqLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<BoussinesqLaw> clone() const;

        //- Selector from dictionary
        inline static autoPtr<BoussinesqLaw> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("Boussinesq");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override
        {
            return true;
        }

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override
        {
            if (T_->isConst())
            {
                return true;
            }
            return false;
        }

        //- Functions for density [kg/m^3]
        defineScalarFieldFunctionsOne(rho, T_);

        //- Functions for compressibility rho/p [s^2/m^2]
        defineScalarFieldFunctionsNone(psi);

        //- Functions for compression factor []
        defineScalarFieldFunctionsNone(Z);

        //- Functions for (Cp - Cv) [J/(kg K)]
        defineScalarFieldFunctionsZero(CpMCv);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update dependent model table
        virtual void updateTable(const word&) override;

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

#include "BoussinesqLawI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
