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
    (c) 2022 Esi Ltd.

Class
    Foam::rhoFunction1

Description
    Incompressible, function1 equation of state.

    Note:
    Compressibility is zero.

SourceFiles
    rhoFunction1I.H
    rhoFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef rhoFunction1_H
#define rhoFunction1_H

#include "memory/autoPtr/autoPtr.H"
#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"
#include "primitives/functions/Polynomial/polynomialFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class rhoFunction1 Declaration
\*---------------------------------------------------------------------------*/

class rhoFunction1
:
    public materialModel,
    public rhoModel,
    public psiModel,
    public ZModel,
    public CpMCvModel
{
    // Private member data

        //- Calculate value
        autoPtr<Function1<scalar>> rhoFunction1_;

        //- Smallest allowed rho
        scalar rhoMin_;

        //- Largest allowed rho
        scalar rhoMax_;


public:

    // Constructors

        //- Construct from dictionary and object registry
        rhoFunction1
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );

        //- Construct and return a clone
        inline autoPtr<rhoFunction1> clone() const;

        //- Selector from dictionary
        inline static autoPtr<rhoFunction1> New
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Member functions

        //- Return the instantiated type name
        TypeName("rhoFunction1");

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
        defineScalarFieldFunctionsTwo(CpMCv, p_, T_);

        //- Return model dependence
        virtual const depList& dep() const override
        {
            return materialModel::dep();
        }

        //- Update dependent model table
        virtual void updateTable(const word&) override
        {}

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

#include "rhoFunction1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
