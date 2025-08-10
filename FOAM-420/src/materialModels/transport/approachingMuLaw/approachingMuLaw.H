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
    (c) 2011 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

Class
    Foam::approachingMuLaw

Description
    A material model which is approaching the desired viscosity mu_desired,
    starting from mu_start.

    A blended algorithm is:
    if (current_time < some_fixed_time_)
    {
        scalar alpha = current_time / some_fixed_time;
        return (alpha * mu_desired_ + (1 - alpha) * mu_start_);
    }
    else
    {
        return mu_desired_;
    }

    However this one behaves better in my tests:
    if (current_time < some_fixed_time_)
    {
        scalar alpha = (mu_start * 1000 / mu_desired) ^ (1/some_fixed_time)
        return (mu_start / (alpha ^ current_time)) + mu_desired;
    }
    else
    {
        return mu_desired_;
    }


SourceFiles
    approachingMuLaw.C

\*---------------------------------------------------------------------------*/

#ifndef approachingMuLaw_H
#define approachingMuLaw_H

#include "materialModels/materialModel/materialModel.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class approachingMuLaw Declaration
\*---------------------------------------------------------------------------*/

class approachingMuLaw
:
    public materialModel,
    public muModel
{
    // Private data

        //- Approaching coefficients
        scalar mu_desired_;
        scalar mu_start_;
        scalar some_fixed_time_;


public:

    //- Runtime type information
    TypeName("approachingNu");


    // Constructors

        //- Construct from dictionary
        approachingMuLaw
        (
            const objectRegistry& obr,
            const dictionary& dict,
            const word& phaseName,
            const word& specieName,
            const word& name
        );


    // Destructor

        ~approachingMuLaw()
        {}


    // Member Functions

        //- Dynamic viscosity [kg/ms]
        defineScalarFieldFunctionsZero(mu);

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

        //- Re-read data
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
