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
    (c) 2021-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fieldLimiter.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldLimiter, 0);
    addToRunTimeSelectionTable(materialModel, fieldLimiter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldLimiter::fieldLimiter
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
    fieldLimiter::read();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::fieldLimiter::castScalarModel
(
    const word& modelName
)
{
    if (modelName == limitModel::typeName)
    {
        return dynamic_cast<limitModel*>(this);
    }
    return nullptr;
}


bool Foam::fieldLimiter::read()
{
    fieldName_ = dict_->lookupOrDefault<word>("fieldName", "T");
    field_ = constructOrReturnRefFieldPtr<scalar>(fieldName_);
    fieldMin_ = dict_->lookup<scalar>("min");
    fieldMax_ = dict_->lookup<scalar>("max");

    return true;
}


// ************************************************************************* //
