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
    (c) 2017-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "approachingMuLaw.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(approachingMuLaw, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        approachingMuLaw,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::approachingMuLaw::approachingMuLaw
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
:
    materialModel(obr, dict, phaseName, specieName, name)
{
    approachingMuLaw::read();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::approachingMuLaw::castScalarModel
(
    const word& modelName
)
{
    if (modelName == muModel::typeName)
    {
        return dynamic_cast<muModel*>(this);
    }
    return nullptr;
}


Foam::scalar Foam::approachingMuLaw::mu() const
{
    const scalar current_time((mesh_.time()).value());

    if (current_time < some_fixed_time_)
    {
        const scalar alpha
            = pow((mu_start_*1000/mu_desired_), (1.0/some_fixed_time_));

        return
        (
            (mu_start_ /(pow(alpha, current_time))) + mu_desired_
        );
    }
    else
    {
        return mu_desired_;
    }
}


bool Foam::approachingMuLaw::read()
{
    mu_desired_ = dict_->lookup<scalar>("mu_desired");
    mu_start_ = dict_->lookup<scalar>("mu_start");
    some_fixed_time_ = dict_->lookup<scalar>("some_fixed_time");

    return true;
}


// ************************************************************************* //
