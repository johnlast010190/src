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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "constSolid.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constSolid, 0);
    addToRunTimeSelectionTable
    (
        materialModel,
        constSolid,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constSolid::constSolid
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
    constSolid::read();
}


Foam::autoPtr<Foam::constSolid>
Foam::constSolid::clone() const
{
    return autoPtr<constSolid>
    (
        new constSolid(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::baseModels<Foam::scalar>* Foam::constSolid::castScalarModel
(
    const word& modelName
)
{
    if (modelName == kappaModel::typeName)
    {
        return dynamic_cast<kappaModel*>(this);
    }
    return nullptr;
}


Foam::baseModels<Foam::vector>* Foam::constSolid::castVectorModel
(
    const word& modelName
)
{
    if (modelName == vKappaModel::typeName)
    {
        return dynamic_cast<vKappaModel*>(this);
    }
    return nullptr;
}


bool Foam::constSolid::read()
{
    isAniso_ = dict_->lookupOrDefault<bool>("anisotropic", false);
    if (dict_->found("kappa") && isAniso_)
    {
        kappa_ = dict_->lookup<vector>("kappa");
        magKappa_ = mag(kappa_);
    }
    else if (dict_->found("kappa"))
    {
        const scalar kappa(dict_->lookup<scalar>("kappa"));
        kappa_ = vector(kappa, kappa, kappa);
        magKappa_ = kappa;
    }
    else
    {
        FatalErrorInFunction
            << "\"kappa\" not found."
            << nl << exit(FatalError);
    }

    return true;
}


// ************************************************************************* //
