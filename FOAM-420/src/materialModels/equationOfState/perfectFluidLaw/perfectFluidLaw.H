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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021-2022 Esi Ltd.

Class
    Foam::perfectFluidLaw

Description
    Perfect fluid equation of state.

SourceFiles
    perfectFluidLawI.H
    perfectFluidLaw.C

\*---------------------------------------------------------------------------*/

#ifndef perfectFluidLaw_H
#define perfectFluidLaw_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class perfectFluidLaw Declaration
\*---------------------------------------------------------------------------*/

class perfectFluidLaw
:
    public materialModel,
    public rhoModel,
    public SDepartureModel,
    public psiModel,
    public ZModel,
    public CpMCvModel
{
    // Private member data

        //- Fluid constant
        scalar R_;

        //- The reference density
        scalar rho0_;

        //- Quick model access
        enum {R, modelsEnumSize_};


public:


    // Constructors

        //- Construct from dictionary and object registry
        perfectFluidLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<perfectFluidLaw> clone() const;

        //- Selector from dictionary
        inline static autoPtr<perfectFluidLaw> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("perfectFluid");

        //- Is the equation of state is incompressible i.e. rho != f(p)
        virtual bool incompressible() const override;

        //- Is the equation of state is isochoric i.e. rho = const
        virtual bool isochoric() const override;

        //- Functions for density [kg/m^3]
        defineScalarFieldFunctionsTwo(rho, p_, T_);

        //- Functions for entropy [J/(kg K)]
        defineScalarFieldFunctionsOne(SDeparture, p_);

        //- Functions for compressibility rho/p [s^2/m^2]
        defineScalarFieldFunctionsOne(psi, T_);

        //- Functions for compression factor []
        defineScalarFieldFunctionsZero(Z);

        //- Functions for (Cp - Cv) [J/(kg K)]
        defineScalarFieldFunctionsNone(CpMCv);

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

#include "perfectFluidLawI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
